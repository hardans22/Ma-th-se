include("__functions__.jl")
include("__printing__.jl")
include("__sub_problem__.jl")

function benders_decomp(env, instance_data, sp_method, output_file, nbr_thread, silent, time_limit)
    #= if graph_reduc
        instance_data = graph_reduction(instance_data)
    end
    =#
    graph       = instance_data.graph
    node_sets   = graph.node_sets
    arc_sets    = graph.arc_sets
    fl_data     = instance_data.fl_data
    other_data  = instance_data.other_data

    

    H_aircraft  = node_sets.H_aircraft
    fl_nodes    = node_sets.fl_nodes
    a_nodes     = node_sets.ac_nodes
    V_wt_st     = vcat(a_nodes, fl_nodes) 

    A   = graph.arcs
    A_S = arc_sets.arcs_S
    A_T = arc_sets.arcs_T
    A_M = arc_sets.arcs_M
    L_M = node_sets.mtn_nodes
    V   = graph.nodes
    d   = fl_data.d

    nbr_K   = instance_data.nbr_ac
    A_M_bar = arc_sets.arcs_M_bar
    
    #= if graph_reduc
        F_bar = other_data["F_bar"]
        d_bar = other_data["d_bar"]
    end
 =#
    pred_A, succ_A = Dict(), Dict()
    pred_A, succ_A = update_neighborhood(A, pred_A, succ_A)
    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(A_M, pred_A_M, succ_A_M)
    pred_A_M_bar, succ_A_M_bar = Dict(), Dict()
    pred_A_M_bar, succ_A_M_bar = update_neighborhood(A_M_bar, pred_A_M_bar, succ_A_M_bar)
    
    # ===================== Master problem =====================
    
    master_model = Model(() -> Gurobi.Optimizer(env))
    GRBsetintparam(env, "OutputFlag", silent ? 0 : 1)
    set_optimizer_attribute(master_model, "Threads", nbr_thread)
    set_optimizer_attribute(master_model, "TimeLimit", time_limit) 
    set_optimizer_attribute(master_model, "LazyConstraints", 1)
    set_optimizer_attribute(master_model, "IntFeasTol", 1e-9)
    set_optimizer_attribute(master_model, "FeasibilityTol", 1e-09)
    set_optimizer_attribute(master_model, "LogFile", output_file)
    set_optimizer_attribute(master_model, "Cuts", 0)
    # Variables
    @variable(master_model, x[arc in A], Bin)
    @variable(master_model, Z[k in a_nodes] >= 0)

    # Contraintes
    #= @constraint(master_model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(master_model, sum(x[i] for i in A_T) == nbr_K)
     =#
    @constraint(master_model, mp_c1[i in V_wt_st], sum(x[(i,j)] for j in get(succ_A, i, [])) == 1)
    @constraint(master_model, mp_c2[i in V_wt_st], sum(x[(j,i)] for j in get(pred_A, i, [])) == sum(x[(i,j)] for j in get(succ_A, i, [])))
    
    @objective(master_model, Min, sum(Z[k] for k in a_nodes))

    
    A_theta2 = [(i,j) for (i,j) in A if j != "t"]
    A_theta3 = [(i,j) for (i,j) in A_M_bar if j != "t"]
    eq_nodes = union(["s"], a_nodes)
    other_data["A_theta2"] = A_theta2
    other_data["A_theta3"] = A_theta3
    other_data["eq_nodes"] = eq_nodes
    other_data["pred_A_M"] = pred_A_M
    
    
    rho_val = Dict(j => 0.0 for j in L_M)
    u_val = Dict(j => 0.0 for j in V)
    y_val = Dict(j => 0.0 for j in L_M) 
         
    nbcp_fais, nbcp_opti, nbr_iter = 0, 0, 0                     
    all_sptime = 0
    best_obj = 0
    best_x = Dict{Tuple{String, String}, Int}
    best_y = Dict{String, Int}
    best_u = Dict{String, Int}
    best_rho = Dict{String, Int}
                     
    sp_model = build_sp(env, instance_data)
    benders_cuts = Vector()
    temp = true
    function my_callback(cb_data)
        status = callback_node_status(cb_data, master_model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end

        nbr_iter += 1
        
        x_val = Dict((i,j) => round(callback_value(cb_data, x[(i,j)])) for (i,j) in A)
        #@assert all(v -> isapprox(v, 0.0, atol=1e-5) || isapprox(v, 1.0, atol=1e-5), values(x_val))
        xx_val = Dict{VariableRef, Float64}(x[a] => round(callback_value(cb_data, x[a])) for a in A)
        #ZZ_val = Dict{VariableRef, Float64}(Z[k] => callback_value(cb_data, Z[k]) for k in a_nodes)
        #xx_val = Dict{VariableRef, Float64}(var => callback_value(cb_data, var) for var in all_variables(master_model))
         
        Z_val = Dict(k => callback_value(cb_data, Z[k]) for k in a_nodes)
        #Z_val = callback_value(cb_data, Z)
        
        aircraft_paths = build_all_paths(x_val, A, a_nodes)

        
        if sp_method == "MILP" 
            result = solve_sp(sp_model, x_val, L_M, V)
        elseif sp_method == "PD"
            result = solve_sp_pd(aircraft_paths, instance_data, a_nodes)
        end 

        sp_obj = result.obj
    
        sub_paths = find_subpath_mtn(result, instance_data, aircraft_paths)
         
        all_sptime += result.time

        if result.status == "OPTIMAL"
            #println("\n\n SP STATUS: ", result.status)
            #=println("Z_VAL = ", Z_val)
            println("SP OBJ = ", sp_obj)
            println("SP TIME = ", result.time)
            print_x_y(x_val, y_val, rho_val, A, L_M, d, H_aircraft)
             =#
            obj_Z = sum(Z_val[k] for k in H_aircraft)
            #= println()
            println(Z_val)
            println("Obj_Z = $obj_Z")
 =#
            gap = sp_obj - obj_Z
            rho_val = result.rho
            u_val = result.u
            y_val = result.y
            
            if gap > 1e-6
                #ADD OPTIMALITY CUT
                #check_cuts(benders_cuts, xx_val; tol=1e-6)

                cut_list = build_opti_cut(sub_paths, x, Z, rho_val)
                for cut in cut_list
                    nbcp_opti += 1
                    if temp
                        push!(benders_cuts, cut)
                    end
                    MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
                end 
                #= cut = build_opti_cut_agg(sub_paths, x, Z, sp_obj)
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
                nbcp_opti += 1
                benders_cuts[nbr_iter] = [cut]
                 =#
                #println(cut_list[1])
                #println("COUPE AJOUTÉE | Z_VAL = ", obj_Z, " | SP = ", sp_obj, " | nb_coupes = ", nbcp_opti)

            else 
                #OPTIMAL SOLUTION FOUND
                best_obj = sp_obj
                best_x = x_val
                best_y = y_val 
                best_u = u_val
                best_rho = rho_val 
                #println() 
                println("\nZ_VAL = ", obj_Z, " | OBJ DU SP = ", result.obj, " | GAP = ", gap)
                #print_x_y(x_val, y_val, rho_val, A, L_M, d, H_aircraft)
                #= if obj_Z == 1579.9999999999975
                    #print_x_y(x_val, y_val, rho_val, A, L_M, d, H_aircraft)
                    println(Z_val)
                    #check_cuts(benders_cuts, xx_val; tol=1e-6)
                end =#
                if gap < -0.1 && temp
                    println("C'est arrivé une fois")
                    temp = false

                    #= println(Z_val)
                    println("Obj_Z = $obj_Z")

                    println("AAAAAALERRRRRRTE IL Y A LE FEU ")
                    #check_cuts(benders_cuts, xx_val; tol=1e-6)
                    print_x_y(x_val, y_val, rho_val, A, L_M, d, a_nodes)
                     =#
                    #= open("Opti_cuts.txt", "w") do fichier
                        for cut in benders_cuts
                            println(fichier, cut)
                        end
                    end =#
                end  
                #= if gap < 0
                    check_cuts(benders_cuts, xx_val; tol=1e-6)
                    print_x_y(x_val, y_val, rho_val, A, L_M, d, H_aircraft)
                end =#  
            end
             
        elseif result.status == "INFEASIBLE" 
            #ADD FEASIBILITY_CUTS
            #println("\nSP STATUS: ", result.status)
            #print_x_y(x_val, y_val, rho_val, A, L_M, d, H_aircraft)
            #l_inf = infeasibility(sp_model)
            #println("\n")
            #println(length(l_inf))
            cut_list = build_feas_cut(sub_paths, x)
            for cut in cut_list
                nbcp_fais += 1  
                if temp
                    push!(benders_cuts, cut)
                end
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
            end 
        else
            println("Erreur: Statut du sous-problème = $(result.status)")
        end
    end
    #write_to_file(master_model, "mp.lp")

    # Set callback
    set_attribute(master_model, MOI.LazyConstraintCallback(), my_callback)
    
    # solve
    optimize!(master_model)
    println("BEST SP OBJ = ", best_obj)

    open("Benders_cuts.txt", "w") do fichier
        for cut in benders_cuts
            println(fichier, cut)
        end
    end

    status = termination_status(master_model)
    nb_coupes = nbcp_fais + nbcp_opti

    #print_x_y(best_x, best_y, best_rho, A, L_M, d, a_nodes)
    obj_val = best_obj
    time = round(solve_time(master_model), digits = 4)
    gap = round(relative_gap(master_model)*100, digits = 2)
    nbr_nodes =  MOI.get(master_model, MOI.NodeCount())
    dual_obj = round(objective_bound(master_model); digits = 2)
    time = round(solve_time(master_model), digits = 4)
    all_sptime = round(all_sptime, digits = 4)
    mp_time = round(time - all_sptime, digits = 4)
    #nbr_mtn = round(Int, sum(y_val))
    mtn_stations_used = []
    
    nbr_mtn = round(Int, sum(values(best_y)))
    
    nbr_sts_used = round(Int, length(mtn_stations_used))
    
    other_info = Dict("gap" => gap, "nbr_nodes" => nbr_nodes, "time" => time, "mp_time" => mp_time,
                    "sp_time" => all_sptime, "dual_obj" => dual_obj, "status" => status, 
                    "nbr_iter" => nbr_iter, "nbr_cuts" => nb_coupes, "nbr_feas_cuts" => nbcp_fais, 
                    "nbr_opti_cuts" => nbcp_opti, "nbr_mtn" => nbr_mtn, "nbr_sts_used" => nbr_sts_used)

    solution = Solution_FH(best_obj, best_x, best_y, best_u, best_rho, other_info)
    
    #= output_results = Dict( "instance" => instance, "obj" => obj_val, "x" => best_x, "y" => best_y, 
                    "u" => best_u, "rho" => best_rho, "gap" => gap, "nbr_nodes" => nbr_nodes, "time" => time, "mp_time" => mp_time,
                    "sp_time" => all_sptime, "dual_obj" => dual_obj,  
                    "status" => status, "nbr_iterations" => nbr_iter, 
                    "nbr_cuts" => nb_coupes, "nbr_feasibility_cuts" => nbcp_fais, 
                    "nbr_optimality_cuts" => nbcp_opti)
 =#
    return solution
end