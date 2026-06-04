include("../../structures.jl")


function model_amrp(env, instance_data, output_file, nbr_thread, silent, graph_reduc, time_limit, work_limit)
    if graph_reduc
        instance_data = graph_reduction(instance_data)
    end

    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data
    other_data = instance_data.other_data

    nbr_K = instance_data.nbr_ac
    max_flt = instance_data.max_flying_time
    
    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes
    L_M = node_sets.mtn_nodes
    L_MS = node_sets.mtn_nodes_stn
    H_aircraft = node_sets.H_aircraft
    V = graph.nodes
    V_wt_st = vcat(a_nodes, fl_nodes)  #Nodes without st nodes

    A_S = arc_sets.arcs_S
    A_T = arc_sets.arcs_T
    A_M = arc_sets.arcs_M
    A_M_bar = arc_sets.arcs_M_bar
    A = graph.arcs

    f = fl_data.init_flying_time
    d = fl_data.d

    MS = other_data["mtn_stations"]
    
    if graph_reduc
        F_bar = other_data["F_bar"]
        d_bar = other_data["d_bar"]
    end

    predecessors, successors = Dict(), Dict()
    predecessors, successors = update_neighborhood(A, predecessors, successors)
    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(A, pred_A_M, succ_A_M)
    
    # ===================== Model =====================
    model = Model(() -> Gurobi.Optimizer(env))
    #GRBsetintparam(env, "OutputFlag", 0)
    if silent 
        set_optimizer_attribute(model, "LogToConsole", 0)
    end 
    set_optimizer_attribute(model, "Threads", nbr_thread)
    set_optimizer_attribute(model, "TimeLimit", time_limit) 
    set_optimizer_attribute(model, "LogFile", output_file)
    #set_optimizer_attribute(model, "WorkLimit", work_limit)
    #set_optimizer_attribute(model, "Presolve", 0)
    #set_optimizer_attribute(model, "Cuts", 2)           # Aggressive

    # ===================== Decision variables =====================
    @variable(model, x[k in A], Bin)                          #Arc (i,j) selected
    @variable(model, y[j in L_M], Bin)                        #1 if a maintenance is performed at node j 

    # ===================== Route Constraints =====================
    #@constraint(model, sum(x[i] for i in A_S) == nbr_K)
    #@constraint(model, sum(x[i] for i in A_T) == nbr_K)
    
    #@constraint(model, sum(x[i] for i in A_T) == sum(x[i] for i in A_S))
    @constraint(model, c1[i in V_wt_st], sum(x[(i,j)] for j in get(successors, i, [])) == 1)
    @constraint(model, c2[i in V_wt_st], sum(x[(j,i)] for j in get(predecessors, i, [])) == sum(x[(i,j)] for j in get(successors, i, [])))
    @constraint(model, c3[j in L_M], y[j] <= sum(x[(i, j)] for i in get(pred_A_M, j, [])))
    #@constraint(model, sum(y[j] for j in L_M) <= length(H_aircraft))
    
    # =========== Variables ===========                        
    if graph_reduc 
        @variable(model, d_bar[j] <= u[j in V] <= F_bar[j])                  #Accumulatad flying time at node j 
    else
        @variable(model, d[j] <= u[j in V] <= max_flt)                  #Accumulatad flying time at node j 
    end 
    @variable(model, rho[j in L_M] >= 0)                              #Remaining flying time at node j

    # ========  Flying time constraints ========
    @constraint(model, c4[(i,j) in A_M], rho[j] >= max_flt*x[(i, j)] - u[i] - (max_flt - d[i])*(1 - y[j]))
    @constraint(model, [(i,j) in A; j != "t"], u[j] <= u[i] + d[j] + (max_flt - d[i] - d[j])*(1 - x[(i,j)]), base_name = "c5[($i,$j)]")
    @constraint(model, c6[j in L_M], u[j] <= max_flt -(max_flt - d[j])*y[j])
    @constraint(model, [(i,j) in A_M_bar; j != "t"], u[j] >= u[i] + d[j] - max_flt*(1 - x[(i,j)]), base_name = "c7[($i,$j)]")
    @constraint(model, c8[(i,j) in A_M], u[j] >= u[i] + d[j] - max_flt*(1 - x[(i,j)]) - max_flt*y[j])
    @constraint(model, u["s"] == 0)
    @constraint(model, c9[k in a_nodes], u[k] == f[k]) 

    # =========== Objective ===========
    @objective(model, Min, sum(rho[j] for j in L_M))

    optimize!(model)

    # ===================== Résultats =====================

    # Contrôle du statut
    status = termination_status(model)
    println("STATUT : ", status)
    if status == MOI.OPTIMAL || status == MOI.FEASIBLE_POINT || status == MOI.TIME_LIMIT || status == MOI.OTHER_LIMIT || status == MOI.INTERRUPTED
        # ===================== Variables =====================
        sx, sy = Dict((i,j) => JuMP.value(x[(i,j)]) for (i,j) in A), Dict(j => JuMP.value(y[j]) for j in L_M)
        su, s_rho = Dict(j => JuMP.value(u[j]) for j in V), Dict(j => round(Int, JuMP.value(rho[j])) for j in L_M)
        
        active_j = Set(j for j in L_M if value(y[j]) >= 0.9)
        obj_rho = isempty(active_j) ? 0 : sum(s_rho[j] for j in active_j)

        obj_val = objective_value(model)
        gap = round(relative_gap(model)*100, digits = 2)
        nbr_nodes =  MOI.get(model, MOI.NodeCount())
        dual_obj = round(objective_bound(model); digits = 2)
        time = round(solve_time(model), digits = 2)
        work_unit = round(MOI.get(model, Gurobi.ModelAttribute("Work")), digits = 2)


        nbr_mtn = round(Int, sum(values(sy)))

        mtn_stations_used = []
        active_j = Set(j for j in L_M if value(y[j]) >= 0.9)
        remaining_stations = setdiff(Set(MS), mtn_stations_used)

        for ms in remaining_stations
            if !isempty(intersect(Set(L_MS[ms]), active_j))
                push!(mtn_stations_used, ms)
            end
        end
        
        nbr_sts_used = round(Int, length(mtn_stations_used))

        other_info = Dict("obj_rho" => obj_rho, "gap" => gap, "nbr_nodes" => nbr_nodes, "status" => status, 
                    "time" => time, "work" => work_unit, "dual_obj" => dual_obj, "nbr_mtn" => nbr_mtn, 
                    "nbr_sts_used" => nbr_sts_used, "mtn_stations_used" => mtn_stations_used)
        
        solution = Solution_FH(obj_val, sx, sy, su, s_rho, other_info)

        return solution
    else
        temp = Dict(j => 0 for j in V)
        other_info = Dict("status" => status)
        return Dict("status" => status)
    end
end
