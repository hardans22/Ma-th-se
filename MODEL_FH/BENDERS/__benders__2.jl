include("../../structures.jl")
include("__cuts_func__.jl")
include("__printing__.jl")
include("__sub_problem__.jl")


function benders_decomp_2(env, instance_data, sp_method, cut_type, output_file, nbr_thread, silent, time_limit)
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
    set_optimizer_attribute(master_model, "IntFeasTol", 1e-9)
    set_optimizer_attribute(master_model, "FeasibilityTol", 1e-09)
    set_optimizer_attribute(master_model, "LogFile", output_file)
    #set_optimizer_attribute(master_model, "Cuts", 0)
    set_optimizer_attribute(master_model, "Symmetry", 2)

    # ===================== Decision Variables =====================

    @variable(master_model, x[arc in A], Bin)

    if cut_type == "desagg"
        @variable(master_model, Z[k in a_nodes] >= 0)
        @objective(master_model, Min, sum(Z[k] for k in a_nodes))
    elseif cut_type == "agg"
        @variable(master_model, Z >= 0)
        @objective(master_model, Min, Z)
    end

    # ===================== Constraints =====================

    @constraint(master_model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(master_model, sum(x[i] for i in A_T) == nbr_K)
    @constraint(master_model, mp_c1[i in V_wt_st], sum(x[(i,j)] for j in get(succ_A, i, [])) == 1)
    @constraint(master_model, mp_c2[i in V_wt_st], sum(x[(j,i)] for j in get(pred_A, i, [])) == sum(x[(i,j)] for j in get(succ_A, i, [])))

    # ===================== Initialization =====================

    rho_val = Dict(j => 0.0 for j in L_M)
    u_val   = Dict(j => 0.0 for j in V)
    y_val   = Dict(j => 0.0 for j in L_M)
    nbr_fais, nbr_opti, nbr_fix, nbr_iter = 0, 0, 0, 0
    all_sptime = 0.0
    best_obj = Inf
    best_x, best_y, best_u, best_rho = Dict{Tuple{String,String}, Float64}(), Dict{String,Float64}(), Dict{String,Float64}(), Dict{String,Float64}()

    sp_model  = build_sp(env, instance_data)
    fix_dict  = Dict(k => 0 for k in H_aircraft)
    added_cuts = Set()

    println()
    start_time = time()
    converged  = false

    # ===================== Benders loop =====================
    iteration = 0
    while !converged
        elapsed = time() - start_time
        if elapsed >= time_limit
            println("Time limit reached before convergence.")
            break
        end
        iteration += 1

        # --- Solve master ---
        set_optimizer_attribute(master_model, "TimeLimit", time_limit - elapsed)
        optimize!(master_model)

        mp_status = termination_status(master_model)

        if mp_status == MOI.INFEASIBLE || mp_status == MOI.INFEASIBLE_OR_UNBOUNDED
            println("Master problem infeasible. Stopping.")
            break
        end

        if mp_status != MOI.OPTIMAL && mp_status != MOI.FEASIBLE_POINT
            println("Master problem ended with status: $mp_status. Stopping.")
            break
        end

        nbr_iter += 1

        x_val = Dict((i,j) => round(value(x[(i,j)])) for (i,j) in A)

        if cut_type == "desagg"
            Z_val  = Dict(k => round(value(Z[k])) for k in a_nodes)
            obj_Z  = sum(Z_val[k] for k in H_aircraft)
        elseif cut_type == "agg"
            Z_val  = round(value(Z))
            obj_Z  = Z_val
        end

        LB = objective_bound(master_model)

        # --- Solve subproblem ---
        aircraft_paths = build_all_paths(x_val, A, a_nodes)

        if sp_method == "MILP"
            result = solve_sp(sp_model, x_val, L_M, V)
        elseif sp_method == "PD"
            result = solve_sp_pd(aircraft_paths, instance_data, a_nodes)
        end

        sp_obj     = result.obj
        sub_paths  = find_subpath_mtn(result, instance_data, aircraft_paths)
        all_sptime += result.time
        println("Iteration : $iteration  || LB : $obj_Z  || sp_obj : $sp_obj")
        if result.status == "OPTIMAL"
            gap = abs(sp_obj - obj_Z)
            rho_val = result.rho
            u_val   = result.u
            y_val   = result.y

            if gap > 1e-6
                # --- Add optimality cut ---
                if cut_type == "desagg"
                    cut_dict = build_opti_cut(sub_paths, x, Z, rho_val, fix_dict, added_cuts)
                    for cut in cut_dict["opti"]
                        JuMP.add_constraint(master_model, cut)
                        nbr_opti += 1
                    end
                    #= for cut in cut_dict["fix"]
                        JuMP.add_constraint(master_model, cut)
                        nbr_fix += 1
                    end =#
                elseif cut_type == "agg"
                    for cut in build_opti_cut_agg(sub_paths, x, Z, sp_obj, rho_val)
                        JuMP.add_constraint(master_model, cut)
                        nbr_opti += 1
                    end
                end
            else
                
                best_obj  = obj_Z
                best_x    = copy(x_val)
                best_y    = copy(y_val)
                best_u    = copy(u_val)
                best_rho  = copy(rho_val)
                println("Converged: UB = $best_obj, LB = $LB")
                converged = true
            end

        elseif result.status == "INFEASIBLE"
            # --- Add feasibility cut ---
            for cut in build_feas_cut(sub_paths, x, d, A_M)
                JuMP.add_constraint(master_model, cut)
                nbr_fais += 1
            end
        else
            println("Erreur: Statut du sous-problème = $(result.status)")
            break
        end
    end

    # ===================== End Benders =====================

    nbr_cuts = nbr_fais + nbr_opti + nbr_fix
    all_time  = round(time() - start_time, digits=4)
    all_sptime = round(all_sptime, digits=4)
    mp_time   = round(all_time - all_sptime, digits=4)
    gap       = best_obj == Inf ? Inf : round((best_obj - objective_bound(master_model)) / (abs(best_obj) + 1e-10) * 100, digits=2)
    nbr_nodes = MOI.get(master_model, MOI.NodeCount())
    dual_obj  = round(objective_bound(master_model), digits=2)

    if best_obj == Inf
        println("================ NO SOLUTION FOUND ================")
        other_info = Dict("gap" => gap, "nbr_nodes" => nbr_nodes, "time" => all_time, "mp_time" => mp_time,
                          "sp_time" => all_sptime, "dual_obj" => dual_obj, "status" => "NO SOLUTION FOUND",
                          "nbr_iter" => nbr_iter, "nbr_cuts" => nbr_cuts, "nbr_feas_cuts" => nbr_fais,
                          "nbr_opti_cuts" => nbr_opti, "nbr_fix_cuts" => nbr_fix)
    else
        println("================ SOLUTION FOUND ================")
        status    = termination_status(master_model)
        nbr_mtn   = round(Int, sum(values(best_y)))
        other_info = Dict("gap" => gap, "nbr_nodes" => nbr_nodes, "time" => all_time, "mp_time" => mp_time,
                          "sp_time" => all_sptime, "dual_obj" => dual_obj, "status" => status,
                          "nbr_mtn" => nbr_mtn, "nbr_iter" => nbr_iter, "nbr_cuts" => nbr_cuts,
                          "nbr_feas_cuts" => nbr_fais, "nbr_opti_cuts" => nbr_opti, "nbr_fix_cuts" => nbr_fix,
                          "nbr_sts_used" => 0)
    end

    return Solution_FH(best_obj, best_x, best_y, best_u, best_rho, other_info)
end