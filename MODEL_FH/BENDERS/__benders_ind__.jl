include("../../structures.jl")
include("__cuts_func_ind__.jl")
include("__printing__.jl")
include("__sub_problem_ind__.jl")


function benders_decomp_ind(env, instance_data, sp_method, output_file, nbr_thread, silent, time_limit)
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

    pred_A, succ_A = Dict(), Dict()
    pred_A, succ_A = update_neighborhood(A, pred_A, succ_A)
    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(A_M, pred_A_M, succ_A_M)
    pred_A_M_bar, succ_A_M_bar = Dict(), Dict()
    pred_A_M_bar, succ_A_M_bar = update_neighborhood(A_M_bar, pred_A_M_bar, succ_A_M_bar)

    other_data["pred_A_M"] = pred_A_M

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

    @variable(master_model, Z[k in a_nodes] >= 0)
    @objective(master_model, Min, sum(Z[k] for k in a_nodes))

    # ===================== Constraints =====================
    @constraint(master_model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(master_model, sum(x[i] for i in A_T) == nbr_K)
    @constraint(master_model, mp_c1[i in V_wt_st], sum(x[(i,j)] for j in get(succ_A, i, [])) == 1)
    @constraint(master_model, mp_c2[i in V_wt_st], sum(x[(j,i)] for j in get(pred_A, i, [])) == sum(x[(i,j)] for j in get(succ_A, i, [])))

    # ===================== Initialization =====================

    nbr_fais, nbr_opti, nbr_fix, nbr_iter = 0, 0, 0, 0
    all_sptime = 0
    all_build_t = 0
    best_obj = Inf
    best_x, best_y = Dict{Tuple{String,String}, Float64}(), Dict{String, Float64}()
    best_u, best_rho = Dict{String, Float64}(), Dict{String, Float64}()

    if sp_method == "MILP"
        sp_model_c = build_sp(env, instance_data)
    end

    fix_dict = Dict(k => 0 for k in H_aircraft)
    obj_Z = nothing

    # ===================== Callback function =====================
    println()
    start_time = time()
    #sp_model = build_sp(env, instance_data)

    function my_callback(cb_data, cb_where::Cint)
        if cb_where != GRB_CB_MIPSOL && cb_where != GRB_CB_MIPNODE
            return
        end

        if cb_where == GRB_CB_MIPNODE
            resultP = Ref{Cint}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
            if resultP[] != GRB_OPTIMAL
                return
            end
        end

        status = callback_node_status(cb_data, master_model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end
        nbr_iter += 1

        Gurobi.load_callback_variable_primal(cb_data, cb_where)
        x_val = Dict((i,j) => round(callback_value(cb_data, x[(i,j)])) for (i,j) in A)
        Z_val = Dict(k => round(callback_value(cb_data, Z[k])) for k in a_nodes)
        obj_Z = sum(Z_val[k] for k in a_nodes)

        aircraft_paths = build_all_paths(x_val, A, a_nodes)
        sets_for_SP = set_for_SP(aircraft_paths, instance_data)

        u_val   = Dict{String, Float64}(j => 0.0 for j in V)
        y_val   = Dict{String, Int}(j => 0 for j in L_M)
        rho_val = Dict{String, Float64}(j => 0.0 for j in L_M)

        sp_obj = 0.0
        # Flag global : reste true tant qu'aucune coupe n'a été nécessaire
        all_satisfied = true
        all_optimal = true
        fixing_cuts = Dict()
        println()
        println("NBR ITERATION : ", nbr_iter)
        
        for ac in a_nodes
            set_for_SP_ac = sets_for_SP[ac]
            if sp_method == "MILP"
                build_t = 0.0
                build_t = @elapsed begin
                    sp_model = build_sp_ind(env, ac, set_for_SP_ac, instance_data)
                end
                result = solve_sp_ind(sp_model, x_val, set_for_SP_ac.rdc_L_M, set_for_SP_ac.rdc_V, set_for_SP_ac.rdc_A)
                #result = solve_sp_ind_(sp_model, x_val, A, L_M, set_for_SP_ac.rdc_V, set_for_SP_ac.rdc_A, set_for_SP_ac.rdc_L_M)
            elseif sp_method == "HT"
                build_t = 0.0
                result = solve_sp_H_ind(ac, set_for_SP_ac, aircraft_paths[ac], instance_data)
            end

            irr_path = irreductible_path_ind(aircraft_paths[ac], set_for_SP_ac, result, instance_data)
            sp_obj += result.obj
            all_sptime += (result.time + build_t)
            all_build_t += build_t
            
            println("AIRCRAFT : ", ac)
            println("STATUS : ", result.status)
            if result.status == "OPTIMAL"
                gap = abs(result.obj - Z_val[ac])
                if gap < 0.0
                    println("PPPPPROBLEMMMMMMMMMME")
                    println(gap)
                end 
                if gap > 1e-6
                    # ADD OPTIMALITY CUT
                    cut = build_opti_cut_ind(ac, irr_path, x, Z, result)
                    #= println()
                    println(ac)
                    println("OPTIMALITY CUT :", cut)
                     =#
                    if cut !== nothing
                        MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
                        nbr_opti += 1
                        all_satisfied = false
                    end
                else
                    # Sous-problème satisfait : on conserve ses infos
                    for i in keys(result.u)
                        u_val[i] = result.u[i]
                    end
                    for j in keys(result.y)
                        y_val[j], rho_val[j] = result.y[j], result.rho[j]
                    end
                end

                if ac in H_aircraft && result.obj == 0 && fix_dict[ac] != 1
                    #println("iteration : ", nbr_iter)
                    cut = build_fix_cut(ac, aircraft_paths[ac], irr_path, x, fix_dict)
                    if cut !== nothing 
                        #= println()
                        println("ON FIXE DES CERTAINES VARIABLES")
                        println() =#
                        fixing_cuts[ac] = cut
                        #MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
                        #nbr_fix += 1
                    end
                end

            elseif result.status == "INFEASIBLE"
                cut = build_feas_cut_ind(irr_path, x)
                #= println()
                println(ac)
                println("NO GOOD CUT :", cut)
                 =#
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
                nbr_fais += 1
                all_satisfied = false
                all_optimal = false
            end 
        end
        if all_optimal
            for cut in values(fixing_cuts)
                println()
                println("ON FIXE DES CERTAINES VARIABLES")
                println()
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
                nbr_fix += 1
            end 
        end
        
        # ---- Mise à jour de la meilleure solution complète ----
        if all_satisfied && obj_Z < best_obj - 1e-6
            best_obj = obj_Z
            best_x   = copy(x_val)
            best_y   = copy(y_val)
            best_u   = copy(u_val)
            best_rho = copy(rho_val)
        end
    end

    # Set callback
    set_attribute(master_model, Gurobi.CallbackFunction(), my_callback)

    # solve
    optimize!(master_model)


    println("Temps de construction des modèles : ", all_build_t)

    # ===================== End Benders =====================
    nbr_cuts = nbr_fais + nbr_opti + nbr_fix
    all_time = round(solve_time(master_model), digits = 4)
    gap = round(relative_gap(master_model)*100, digits = 2)
    nbr_nodes = MOI.get(master_model, MOI.NodeCount())
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
        nbr_mtn = round(Int, sum(values(best_y)))
        nbr_sts_used = 0   # à corriger si tu veux vraiment calculer les stations utilisées

        other_info = Dict("gap" => gap, "nbr_nodes" => nbr_nodes, "time" => all_time, "mp_time" => mp_time,
                        "sp_time" => all_sptime, "dual_obj" => dual_obj, "status" => status,
                        "nbr_mtn" => nbr_mtn, "nbr_iter" => nbr_iter, "nbr_cuts" => nbr_cuts,
                        "nbr_feas_cuts" => nbr_fais, "nbr_opti_cuts" => nbr_opti, "nbr_fix_cuts" => nbr_fix,
                        "nbr_sts_used" => nbr_sts_used)
    end

    solution = Solution_FH(best_obj, best_x, best_y, best_u, best_rho, other_info)

    return solution
end