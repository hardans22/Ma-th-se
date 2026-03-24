include("functions_for_benders.jl")
include("__sub_problem__.jl")

function benders_decomp(env, instance_data, output_file, nbr_thread, silent, time_limit)
    if graph_reduc
        instance_data = graph_reduction(instance_data)
    end

    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data
    other_data = instance_data.other_data

    nbr_K = instance_data.nbr_ac
    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes
    L_M = node_sets.mtn_nodes
    V = graph.nodes
    V_wt_st = vcat(a_nodes, fl_nodes) 

    A = graph.arcs
    A_S = arc_sets.arcs_S
    A_T = arc_sets.arcs_T
    A_M = arc_sets.arcs_M
    A_M_bar = arc_sets.arcs_M_bar
    
    if graph_reduc
        F_bar = other_data["F_bar"]
        d_bar = other_data["d_bar"]
    end

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
    set_optimizer_attribute(master_model, "LogFile", output_file)
    #set_optimizer_attribute(master_model, "Heuristics", 0)

    # Variables
    @variable(master_model, x[k in A], Bin)
    @variable(master_model, Z >= 0)

    # Contraintes
    #= @constraint(master_model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(master_model, sum(x[i] for i in A_T) == nbr_K)
     =#
    @constraint(master_model, mp_c1[i in V_wt_st], sum(x[(i,j)] for j in get(succ_A, i, [])) == 1)
    @constraint(master_model, mp_c2[i in V_wt_st], sum(x[(j,i)] for j in get(pred_A, i, [])) == sum(x[(i,j)] for j in get(succ_A, i, [])))
    
    @objective(master_model, Min, Z)

    
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

    nbcp_fais, nbcp_opti, nbcp_cover, nb_iterations = 0, 0, 0, 0
    sptime_accumulated = 0
    last_sp_obj = 0

    sp_model = build_sp(env, instance_data)
    constraints_list = []
    function my_callback(cb_data)
        status = callback_node_status(cb_data, master_model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end

        nb_iterations += 1
        
        x_val = Dict((i,j) => callback_value(cb_data, x[(i,j)]) for (i,j) in A)
        Z_val = callback_value(cb_data, Z)
        
        result = solve_sp(sp_model, x_val, L_M)
        sp_obj = result.obj
        sub_paths = find_subpath_mtn(x_val, result, instance_data)
        
        sptime_accumulated += result.time

        if result.status == MOI.OPTIMAL
            #= println("\n\n SP STATUS: ", result.status)
            println("Z_VAL = ", Z_val)
            println("SP OBJ = ", sp_obj)
            println("SP TIME = ", result.time)
            print_x_y(x_val, y_val, rho_val, A, L_M, d, H_aircraft)
             =#
            y_val = result.y
            rho_val = result.rho
            gap = sp_obj - Z_val
        
            if gap > 1e-6
                #ADD OPTIMALITY CUT
                nbcp_opti += 1
                cut_list = build_opti_cut(instance_data, sub_paths, x, Z, rho_val)
                for cut in cut_list
                    MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut)
                end 
                println(cut_list[1])
                println("COUPE AJOUTÉE | Z_VAL = ", Z_val, " | SP = ", sp_obj, " | nb_coupes = ", nbcp_opti)
            else 
                #OPTIMAL SOLUTION FOUND
                last_sp_obj = sp_obj
                #println() 
                println("Z_VAL = ", Z_val, " | OBJ DU SP = ", result.obj, " | GAP = ", gap)
            end
            
        elseif result.status == MOI.INFEASIBLE 
            #ADD FEASIBILITY_CUTS
            #println("\nSP STATUS: ", result.status)
            #print_x_y(x_val, y_val, rho_val, A, L_M, d, H_aircraft)
            nbcp_fais += 1   
            cut_list = build_feas_cut(sub_paths, x)
            for cut in cut_list
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
    println("LAST SP OBJ = ", last_sp_obj)

    status = termination_status(master_model)
    nb_coupes = nbcp_fais + nbcp_opti + nbcp_cover
    x_val = value.(x)

    #print_x_y(x_val, y_val, rho_val, A, L_M, d, a_nodes)
    obj_val = objective_value(master_model)
    time = round(solve_time(master_model), digits = 2)
    gap = round(relative_gap(master_model)*100, digits = 2)
    nbr_nodes =  MOI.get(master_model, MOI.NodeCount())
    dual_obj = round(objective_bound(master_model); digits = 2)
    time = round(solve_time(master_model), digits = 2)
    sptime_accumulated = round(sptime_accumulated, digits = 2)
    mp_time = round(time - sptime_accumulated, digits = 2)
    #nbr_mtn = round(Int, sum(y_val))
    mtn_stations_used = []
    
    nbr_sts_used = round(Int, length(mtn_stations_used))

    output_results = Dict( "instance" => instance, "obj" => obj_val, "x" => x_val, "y" => y_val, 
                    "u" => u_val, "rho" => rho_val, "gap" => gap, "nbr_nodes" => nbr_nodes, "time" => time, "mp_time" => mp_time,
                    "sp_time" => sptime_accumulated, "dual_obj" => dual_obj,  
                    "status" => status, "nbr_iterations" => nb_iterations, 
                    "nbr_cuts" => nb_coupes, "nbr_feasibility_cuts" => nbcp_fais, 
                    "nbr_optimality_cuts" => nbcp_opti, "nbr_cover_cuts" => nbcp_cover)

    return output_results
end