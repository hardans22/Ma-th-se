include("__functions__.jl")
include("__printing__.jl")
include("__sub_problem__.jl")

function benders_decomp(env, instance_data, sp_method, cut_type, output_file, nbr_thread, silent, time_limit)
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
    max_flt = instance_data.max_flying_time

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
    
    other_data["pred_A_M"] =  pred_A_M

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

    # ===================== Decision Variables =====================
    
    @variable(master_model, x[arc in A], Bin)

    if cut_type == "desagg"
        @variable(master_model, Z[k in a_nodes] >= 0)
        @objective(master_model, Min, sum(Z[k] for k in a_nodes))
    elseif cut_type == "agg"
        @variable(master_model, Z>= 0)
        @objective(master_model, Min, Z)
    end 

    # ===================== Constraints =====================
    @constraint(master_model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(master_model, sum(x[i] for i in A_T) == nbr_K)
    @constraint(master_model, mp_c1[i in V_wt_st], sum(x[(i,j)] for j in get(succ_A, i, [])) == 1)
    @constraint(master_model, mp_c2[i in V_wt_st], sum(x[(j,i)] for j in get(pred_A, i, [])) == sum(x[(i,j)] for j in get(succ_A, i, [])))
    
    # ===================== Initialization =====================

    rho_val, u_val, y_val = Dict(j => 0.0 for j in L_M), Dict(j => 0.0 for j in V), Dict(j => 0.0 for j in L_M) 
    nbr_fais, nbr_opti, nbr_fix, nbr_iter = 0, 0, 0, 0                     
    all_sptime = 0
    best_obj = 0
    best_x, best_y, best_u, best_rho = Dict{Tuple{String, String}, Int}, Dict{String, Int}, Dict{String, Int}, Dict{String, Int}
                     
    sp_model = build_sp(env, instance_data)
    fix_dict = Dict(k => 0 for k in H_aircraft)
    added_cuts = Set()

    # ===================== Callbacks function =====================

    function my_callback(cb_data)
        status = callback_node_status(cb_data, master_model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end

        nbr_iter += 1
        
        x_val = Dict((i,j) => round(callback_value(cb_data, x[(i,j)])) for (i,j) in A)
        
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
            #=println("\n\n SP STATUS: ", result.status)
            println("Z_VAL = ", obj_Z, " | SP = ", sp_obj, ")
            print_x_y(x_val, y_val, rho_val, A, L_M, d, H_aircraft)=#
             
            if cut_type == "desagg"
                Z_val = Dict(k => callback_value(cb_data, Z[k]) for k in a_nodes)
                obj_Z = sum(Z_val[k] for k in H_aircraft)
            elseif cut_type == "agg"
                Z_val = callback_value(cb_data, Z)
                obj_Z = Z_val
            end 

            gap = abs(sp_obj - obj_Z)
            rho_val = result.rho
            u_val = result.u
            y_val = result.y
            
            if gap > 1e-6
                #ADD OPTIMALITY CUT
                if cut_type == "desagg"
                    cut_dict = build_opti_cut(sub_paths, x, Z, rho_val, fix_dict, added_cuts)
                    opti_cuts = cut_dict["opti"]
                    for cut in opti_cuts
                        MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
                        nbr_opti += 1
                    end

                    fix_cuts = cut_dict["fix"]
                    for cut in fix_cuts
                        MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
                        nbr_fix += 1
                    end

                elseif cut_type == "agg"
                    cut_list = build_opti_cut_agg(sub_paths, x, Z, sp_obj, rho_val)
                    for cut in cut_list 
                        MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
                    end 
                end 
                #println("COUPE AJOUTÉE | Z_VAL = ", obj_Z, " | SP = ", sp_obj, " | nbr_cuts = ", nbr_opti)
            else 
                best_obj, best_x, best_y, best_u, best_rho = sp_obj, x_val, y_val, u_val, rho_val
                #println("\nZ_VAL = ", obj_Z, " | OBJ DU SP = ", result.obj, " | GAP = ", gap)
            end
             
        elseif result.status == "INFEASIBLE" 
            #ADD FEASIBILITY_CUTS
            #println("\nSP STATUS: ", result.status)
            cut_list = build_feas_cut(sub_paths, x, d, A_M)
            for cut in cut_list
                nbr_fais += 1  
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

    # ===================== End Benders =====================

    status = termination_status(master_model)
    nbr_cuts = nbr_fais + nbr_opti + nbr_fix

    time = round(solve_time(master_model), digits = 4)
    gap = round(relative_gap(master_model)*100, digits = 2)
    nbr_nodes =  MOI.get(master_model, MOI.NodeCount())
    dual_obj = round(objective_bound(master_model); digits = 2)
    time = round(solve_time(master_model), digits = 4)
    all_sptime = round(all_sptime, digits = 4)
    mp_time = round(time - all_sptime, digits = 4)
    mtn_stations_used = []
    
    nbr_mtn = round(Int, sum(values(best_y)))
    
    nbr_sts_used = round(Int, length(mtn_stations_used))
    
    other_info = Dict("gap" => gap, "nbr_nodes" => nbr_nodes, "time" => time, "mp_time" => mp_time,
                    "sp_time" => all_sptime, "dual_obj" => dual_obj, "status" => status, 
                    "nbr_mtn" => nbr_mtn, "nbr_iter" => nbr_iter, "nbr_cuts" => nbr_cuts, 
                    "nbr_feas_cuts" => nbr_fais, "nbr_opti_cuts" => nbr_opti, "nbr_fix_cuts" => nbr_fix, 
                    "nbr_sts_used" => nbr_sts_used)

    solution = Solution_FH(best_obj, best_x, best_y, best_u, best_rho, other_info)
    
    return solution
end