include("__functions__.jl")
include("__printing__.jl")
include("__sub_problem__.jl")

function benders_no_callback(env, instance_data, sp_method, output_file, nbr_thread, silent, time_limit)

    # ===================== DATA =====================
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
    A_M = arc_sets.arcs_M
    L_M = node_sets.mtn_nodes
    V   = graph.nodes
    A_M_bar = arc_sets.arcs_M_bar


    # ===================== MASTER =====================
    master_model = Model(() -> Gurobi.Optimizer(env))

    set_optimizer_attribute(master_model, "OutputFlag", silent ? 0 : 1)
    set_optimizer_attribute(master_model, "Threads", nbr_thread)
    set_optimizer_attribute(master_model, "TimeLimit", time_limit)
    set_optimizer_attribute(master_model, "Cuts", 0)

    @variable(master_model, x[arc in A], Bin)
    @variable(master_model, Z[k in a_nodes] >= 0)

    pred_A, succ_A = Dict(), Dict()
    pred_A, succ_A = update_neighborhood(A, pred_A, succ_A)
    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(A_M, pred_A_M, succ_A_M)
    pred_A_M_bar, succ_A_M_bar = Dict(), Dict()
    pred_A_M_bar, succ_A_M_bar = update_neighborhood(A_M_bar, pred_A_M_bar, succ_A_M_bar)

    other_data["pred_A_M"] = pred_A_M

    @constraint(master_model, [i in V_wt_st],
        sum(x[(i,j)] for j in get(succ_A, i, [])) == 1)

    @constraint(master_model, [i in V_wt_st],
        sum(x[(j,i)] for j in get(pred_A, i, [])) ==
        sum(x[(i,j)] for j in get(succ_A, i, [])))

    @objective(master_model, Min, sum(Z[k] for k in a_nodes))

    # ===================== INIT =====================
    sp_model = build_sp(env, instance_data)

    max_iter = 1000
    tol = 1e-6

    best_obj = Inf
    best_x = nothing
    best_y = nothing
    best_u = nothing
    best_rho = nothing

    nbcp_fais = 0
    nbcp_opti = 0
    nbr_iter = 0

    all_sptime = 0.0

    # ===================== LOOP =====================
    for it in 1:max_iter

        println("\n--- ITERATION $it ---")

        optimize!(master_model)

        if termination_status(master_model) != MOI.OPTIMAL
            println("Master non optimal")
            break
        end

        nbr_iter += 1

        x_val = Dict((i,j) => round(value(x[(i,j)])) for (i,j) in A)
        Z_val = Dict(k => value(Z[k]) for k in a_nodes)

        aircraft_paths = build_all_paths(x_val, A, a_nodes)

        # ===================== SP =====================
        if sp_method == "MILP"
            result = solve_sp(sp_model, x_val, L_M, V)
        else
            result = solve_sp_pd(aircraft_paths, instance_data, a_nodes)
        end

        all_sptime += result.time
        sp_obj = result.obj
        obj_Z = sum(Z_val[k] for k in H_aircraft)

        println("Master = $obj_Z | SP = $sp_obj")

        # ===================== CAS OPTIMAL =====================
        if result.status == "OPTIMAL"

            gap = sp_obj - obj_Z

            if gap > tol
                println("Ajout coupe d'optimalité")

                sub_paths = find_subpath_mtn(result, instance_data, aircraft_paths)
                rho_val = result.rho

                cut_list = build_opti_cut(sub_paths, x, Z, rho_val)
                for cut in cut_list
                    JuMP.add_constraint(master_model, cut)
                    nbcp_opti += 1
                end

            else
                println("✅ CONVERGENCE")

                best_obj = sp_obj
                best_x = x_val
                best_y = result.y
                best_u = result.u
                best_rho = result.rho

                break
            end

        # ===================== CAS INFEASIBLE =====================
        elseif result.status == "INFEASIBLE"

            println("Ajout coupe de faisabilité")

            sub_paths = find_subpath_mtn(result, instance_data, aircraft_paths)
            cut_list = build_feas_cut(sub_paths, x)

            for cut in cut_list
                JuMP.add_constraint(master_model, cut)
                nbcp_fais += 1
            end

        else
            error("Statut SP inconnu")
        end
    end

    println("\n===== FIN =====")
    println("Obj = $best_obj")

    solution = Solution_FH(
        best_obj,
        best_x,
        best_y,
        best_u,
        best_rho,
        Dict(
            "nbr_iter" => nbr_iter,
            "nbr_cuts" => nbcp_fais + nbcp_opti,
            "nbr_feas_cuts" => nbcp_fais,
            "nbr_opti_cuts" => nbcp_opti,
            "sp_time" => all_sptime
        )
    )

    return solution
end