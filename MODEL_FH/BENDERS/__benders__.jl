include("../../structures.jl")
include("__cuts_func__.jl")
include("__cuts_func_ind__.jl")
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
    #set_attribute(master_model, "MIPGap", 0.0)       # gap relatif → 0
    #set_attribute(master_model, "MIPGapAbs", 0.0)    # gap absolu → 0
    #= set_optimizer_attribute(master_model, "Cutoff", Inf)  # ignore les solutions du master sans validation SP
    set_optimizer_attribute(master_model, "MIPGapAbs", Inf)
     =#
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
    best_obj = Inf
    best_x, best_y, best_u, best_rho = Dict{Tuple{String,String}, Float64}(), Dict{String, Float64}(), Dict{String, Float64}(), Dict{String, Float64}()
                     
    if sp_method == "MILP"
        sp_model = build_sp(env, instance_data)
    end
    fix_dict = Dict(k => 0 for k in H_aircraft)
    added_cuts = Set()
    last_print = 0.0 
    # ===================== Callback function =====================
    println()
    start_time = time() 
    function my_callback(cb_data, cb_where::Cint)
        #= status = callback_node_status(cb_data, master_model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end
         =#

        if cb_where != Gurobi.GRB_CB_MIPSOL
            return
        end 
        nbr_iter += 1

        lb_ref = Ref{Cdouble}(0.0)
        Gurobi.GRBcbget(cb_data, cb_where, Gurobi.GRB_CB_MIPSOL_OBJBND, lb_ref)
        LB = round(lb_ref[], digits = 2)        

        Gurobi.load_callback_variable_primal(cb_data, cb_where)
        x_val = Dict((i,j) => round(callback_value(cb_data, x[(i,j)])) for (i,j) in A)
        if cut_type == "desagg"
            Z_val = Dict(k => round(callback_value(cb_data, Z[k])) for k in a_nodes)
            obj_Z = sum(Z_val[k] for k in H_aircraft)
        elseif cut_type == "agg"
            Z_val = round(callback_value(cb_data, Z))
            obj_Z = Z_val
        end 

        runtime_ref = Ref{Cdouble}(0.0)
        Gurobi.GRBcbget(cb_data, cb_where, Gurobi.GRB_CB_RUNTIME, runtime_ref)
        current_time = runtime_ref[]
        
        #= if current_time - last_print >= 3.0
            last_print  = current_time
            ub_str = best_obj == Inf ? "N/A" : string(best_obj)
            println("INCUMBENT = $ub_str | LB = $LB | OPT_GAP = $(format_gap(best_obj, LB)) | TIME = $(round(current_time))")
        end
        =#
        #= # ================ ENDING CRITERIA ================
        if best_obj - LB <= 1e-4
            println("ENDING")
            #println("INCUMBENT = $ub_str | LB = $LB | OPT_GAP = $(format_gap(best_obj, LB)) | TIME = $(round(current_time))")
            Gurobi.GRBterminate(backend(master_model))
            return
        end 
        =#
        aircraft_paths = build_all_paths(x_val, A, a_nodes) 
        
        if sp_method == "MILP" 
            result = solve_sp(sp_model, x_val, L_M, V)
            #= result_pd = solve_sp_pd(aircraft_paths, instance_data, a_nodes)
            if result_pd.obj != result.obj
                println("DIFFÉRENCE")
                println("Result du milp = ", result.obj)
                println("Result du pd = ", result_pd.obj)
            end =#
        elseif sp_method == "PD"
            result = solve_sp_pd(aircraft_paths, instance_data, a_nodes)
            #= result_milp = solve_sp(sp_model, x_val, L_M, V)
            if result_milp.obj != result.obj
                println("DIFFÉRENCE")
                println("Result du milp = ", result_milp.obj)
                println("Result du pd = ", result.obj)
            end =#
        end 

        sp_obj = result.obj
        sub_paths = find_subpath_mtn(result, instance_data, aircraft_paths)     
        all_sptime += result.time
    
        if result.status == "OPTIMAL"
            #=println("\n\n SP STATUS: ", result.status)
            println("Z_VAL = ", obj_Z, " | SP = ", sp_obj, ")
            print_x_y(x_val, y_val, rho_val, A, L_M, d, H_aircraft)=#  
            
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
                if obj_Z < best_obj
                    best_obj, best_x, best_y, best_u, best_rho = obj_Z, copy(x_val), copy(y_val), copy(u_val), copy(rho_val)
                end

                #println("\nZ_VAL = ", obj_Z, " | OBJ DU SP = ", result.obj, " | GAP_MP_SP = ", gap)
                #println("*INCUMBENT = $best_obj | LB = $LB | OPT_GAP = $(format_gap(best_obj, LB)) | TIME = $(round(current_time))")
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
    #set_attribute(master_model, MOI.LazyConstraintCallback(), my_callback)
    set_attribute(master_model, Gurobi.CallbackFunction(), my_callback)
    # solve
    optimize!(master_model)

    # ===================== End Benders =====================
    nbr_cuts = nbr_fais + nbr_opti + nbr_fix
    all_time = round(solve_time(master_model), digits = 4)
    gap = round(relative_gap(master_model)*100, digits = 2)
    nbr_nodes =  MOI.get(master_model, MOI.NodeCount())
    dual_obj = round(objective_bound(master_model); digits = 2)
    all_sptime = round(all_sptime, digits = 4)
    mp_time = round(all_time - all_sptime, digits = 4)
    
    if best_obj == Inf
        println("================ NO SOLUTION FOUND ================")
        other_info = Dict("gap" => gap, "nbr_nodes" => nbr_nodes, "time" => all_time, "mp_time" => mp_time,
                        "sp_time" => all_sptime, "dual_obj" => dual_obj, "status" => "NO SOLUTION FOUND", 
                        "nbr_iter" => nbr_iter, "nbr_cuts" => nbr_cuts, "nbr_feas_cuts" => nbr_fais, 
                        "nbr_opti_cuts" => nbr_opti, "nbr_fix_cuts" => nbr_fix)
    
    else 
        println("================ SOLUTION FOUND ================")

        status = termination_status(master_model)
        mtn_stations_used = []
        nbr_mtn = round(Int, sum(values(best_y)))
        nbr_sts_used = round(Int, length(mtn_stations_used))
        
        other_info = Dict("gap" => gap, "nbr_nodes" => nbr_nodes, "time" => all_time, "mp_time" => mp_time,
                        "sp_time" => all_sptime, "dual_obj" => dual_obj, "status" => status, 
                        "nbr_mtn" => nbr_mtn, "nbr_iter" => nbr_iter, "nbr_cuts" => nbr_cuts, 
                        "nbr_feas_cuts" => nbr_fais, "nbr_opti_cuts" => nbr_opti, "nbr_fix_cuts" => nbr_fix, 
                        "nbr_sts_used" => nbr_sts_used)
    end

    solution = Solution_FH(best_obj, best_x, best_y, best_u, best_rho, other_info)
    
    return solution
end
