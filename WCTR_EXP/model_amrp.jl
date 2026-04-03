using Dates

include("../structures.jl")

function model_amrp(env, instance_data, output_file, nbr_thread, time_limit, obj_min, silent, graph_reduc)
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
    max_tk = instance_data.max_takeoff
    max_day = instance_data.max_flying_day
    
    st_nodes = node_sets.dummy_nodes
    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes
    L_M = node_sets.mtn_nodes
    L_MS = node_sets.mtn_nodes_stn
    H_aircraft = node_sets.H_aircraft
    NH_aircraft = node_sets.NH_aircraft
    V = graph.nodes
    V_wt_st = vcat(a_nodes, fl_nodes)  #Nodes without st nodes

    A_S = arc_sets.arcs_S
    A_K = arc_sets.arcs_K
    A_F = arc_sets.arcs_F
    A_T = arc_sets.arcs_T
    A_M = arc_sets.arcs_M
    A_M_bar = arc_sets.arcs_M_bar
    A = graph.arcs

    f = fl_data.init_flying_time
    h = fl_data.init_takeoff
    g = fl_data.init_flying_day
    b = fl_data.b
    d = fl_data.d
    tk = fl_data.tk

    fh_day, fh_tk = other_data["fh_day"], other_data["fh_tk"]
    ms_capacity = other_data["ms_capacity"]
    MS = other_data["mtn_stations"]
    M_FL_O, M_FL_D = other_data["M_FL_O"], other_data["M_FL_D"]
    fl_day = other_data["fl_day"]
    nbr_TP = other_data["nbr_TP"]
    TP = 1:nbr_TP
    
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
    set_optimizer_attribute(model, "Threads", nbr_thread)
    set_optimizer_attribute(model, "TimeLimit", time_limit) 
    set_optimizer_attribute(model, "LogFile", output_file)
    #set_optimizer_attribute(model, "LogToConsole", 0)
    #set_optimizer_attribute(model, "Presolve", 0)
    #set_optimizer_attribute(model, "Cuts", 2)           # Aggressive

    if silent
        set_silent(model)
    end 
    
    
    ft_obj, tk_obj, fd_obj = 0, 0, 0

    # ===================== Decision variables =====================
    @variable(model, x[k in A], Bin)                          #Arc (i,j) selected
    @variable(model, y[j in L_M], Bin)                        #1 if a maintenance is performed at node j                        
    if graph_reduc 
        @variable(model, d_bar[j] <= u[j in V] <= F_bar[j])                  #Accumulatad flying time at node j 
    else
        @variable(model, d[j] <= u[j in V] <= max_flt)                  #Accumulatad flying time at node j 
    end 
    @variable(model, rho[j in L_M] >= 0)                              #Remaining flying time at node j
    @variable(model, 0 <= v[j in V] <= max_tk)                      #Accumulatad number of takeoff at node j
    @variable(model, lambda[j in L_M] >= 0, Int)                         #Remaining number of takeoff time at node j
    @variable(model, z[m in MS, t in TP])                     #Number of maintenance opération at station m in period t
    @variable(model, 0 <= w[j in V] <= max_day)                     #Accumulatad number of flying day at node j
    @variable(model, phi[j in L_M] >= 0, Int)                            #Remaining number of flying day at node j
        
    # ===================== Route Constraints =====================
    @constraint(model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(model, sum(x[i] for i in A_T) == nbr_K)
    @constraint(model, c1[i in fl_nodes], sum(x[(i,j)] for j in get(successors, i, [])) == 1)
    @constraint(model, c2[i in V_wt_st], sum(x[(j,i)] for j in get(predecessors, i, [])) == sum(x[(i,j)] for j in get(successors, i, [])))
    @constraint(model, c3[j in L_M], y[j] <= sum(x[(i, j)] for i in get(pred_A_M, j, [])))
    @constraint(model, sum(y[j] for j in L_M) <= length(H_aircraft))

    #Maintenance capacity
    @constraint(model, c24[ms in MS, t in 1:nbr_TP], z[ms,t] == sum(y[i] for i in L_MS[ms] if fl_day[i] == t))
    @constraint(model, c25[ms in MS, t in 1:nbr_TP],  z[ms,t] <= ms_capacity[ms][t]) 
    
    # ========  Flying time constraints ========
    @constraint(model, c4[(i,j) in A_M], rho[j] >= max_flt*x[(i, j)] - u[i] - (max_flt - d[i])*(1 - y[j]))
    @constraint(model, [(i,j) in A; j != "t"], u[j] <= u[i] + d[j] + (max_flt - d[i] - d[j])*(1 - x[(i,j)]), base_name = "c5[($i,$j)]")
    @constraint(model, c6[j in L_M], u[j] <= max_flt -(max_flt - d[j])*y[j])
    @constraint(model, [(i,j) in A_M_bar; j != "t"], u[j] >= u[i] + d[j] - max_flt*(1 - x[(i,j)]), base_name = "c7[($i,$j)]")
    @constraint(model, c8[(i,j) in A_M], u[j] >= u[i] + d[j] - max_flt*(1 - x[(i,j)]) - max_flt*y[j])
    @constraint(model, u["s"] == 0)
    @constraint(model, c9[k in a_nodes], u[k] == f[k]) 
    
    # ========= Takeoff constraints ========
    @constraint(model, c10[(i,j) in A_M], lambda[j] >=  max_tk*x[(i, j)] - v[i] - (max_tk - tk[i])*(1 - y[j]))
    @constraint(model, [(i,j) in A; j != "t"], v[j] <= v[i] + tk[j] + (max_tk - tk[i] - tk[j])*(1 - x[(i,j)]), base_name = "c11[($i,$j)]")
    @constraint(model, c12[j in L_M], v[j] <= max_tk -(max_tk - tk[j])*y[j])
    @constraint(model, [(i,j) in A_M_bar; j != "t"], v[j] >= v[i] + tk[j] - max_tk*(1 - x[(i,j)]), base_name = "c13[($i,$j)]")
    @constraint(model, c14[(i,j) in A_M], v[j] >= v[i] + tk[j] - max_tk*(1 - x[(i,j)]) - max_tk*y[j])
    @constraint(model, v["s"] == 0)
    @constraint(model, c15[k in a_nodes], v[k] == h[k])     

    # ========= Flying day constraints =========
    @constraint(model, c17[(i,j) in A_M], phi[j] >= max_day*x[(i, j)] - w[i] - (max_day - 1)*(1 - y[j]))
    @constraint(model, [(i,j) in A; j != "t"], w[j] <= w[i] + b[(i,j)] + (max_day - b[(i,j)] - 1)*(1 - x[(i,j)]), base_name = "c18[($i,$j)]")
    @constraint(model, c19[j in L_M], w[j] <= max_day - (max_day - 1)*y[j])
    @constraint(model, [(i,j) in A_M_bar; j != "t"], w[j] >= w[i] + b[(i,j)] - max_day*(1 - x[(i,j)]), base_name = "c20[($i,$j)]")
    @constraint(model, c21[(i,j) in A_M], w[j] >= w[i] + b[(i,j)] - max_day*(1 - x[(i,j)]) - max_day*y[j])
    @constraint(model, w["s"] == 0)
    @constraint(model, c22[k in a_nodes], w[k] == g[k])
    @constraint(model, c23[j in V_wt_st], 1 <= w[j])      
    
    # =========== Objective ===========
    if obj_min == "FH"
        @objective(model, Min, sum(rho[j] for j in L_M))
    elseif obj_min == "FC"
        @objective(model, Min, sum(lambda[j] for j in L_M))
    elseif obj_min == "DY"
        @objective(model, Min, sum(phi[j] for j in L_M))
    end
    
    #write_to_file(model, "mon_modele.lp") 
    
    # ===================== Solve =====================
    optimize!(model)

    # ===================== Résultats =====================

    # Contrôle du statut
    status = termination_status(model)
    if status == MOI.OPTIMAL || status == MOI.FEASIBLE_POINT || status == MOI.TIME_LIMIT || status == MOI.INTERRUPTED
        sx, sy = Dict((i,j) => JuMP.value(x[(i,j)]) for (i,j) in A), Dict(j => JuMP.value(y[j]) for j in L_M)
        su, s_rho = Dict(j => JuMP.value(u[j]) for j in V), Dict(j => round(Int, JuMP.value(rho[j])) for j in L_M)
        sv, s_lambda = Dict(j => JuMP.value(v[j]) for j in V), Dict(j => round(Int, JuMP.value(lambda[j])) for j in L_M)
        sw, s_phi = Dict(j => JuMP.value(w[j]) for j in V), Dict(j => round(Int, JuMP.value(phi[j])) for j in L_M)
        sz = JuMP.value.(z)
        
        obj_val = objective_value(model)
        gap = round(relative_gap(model)*100, digits = 2)
        nbr_nodes =  MOI.get(model, MOI.NodeCount())
        dual_obj = round(objective_bound(model); digits = 2)
        time = round(solve_time(model), digits = 2)

        active_j = Set(j for j in L_M if value(y[j]) >= 0.9)
        obj_rho = isempty(active_j) ? 0 : sum(s_rho[j] for j in active_j)
        obj_lambda = isempty(active_j) ? 0 : sum(s_lambda[j] for j in active_j)
        obj_phi = isempty(active_j) ? 0 : sum(s_phi[j] for j in active_j)


        nbr_mtn = round(Int, sum(values(sy)))

        mtn_stations_used = []
        remaining_stations = setdiff(Set(MS), mtn_stations_used)

        for ms in remaining_stations
            if !isempty(intersect(Set(L_MS[ms]), active_j))
                push!(mtn_stations_used, ms)
            end
        end
        
        nbr_sts_used = round(Int, length(mtn_stations_used))

        other_info = Dict("sz" => sz, "obj_rho" => obj_rho,"obj_lambda" => obj_lambda, "obj_phi" => obj_phi, 
                    "gap" => gap, "nbr_nodes" => nbr_nodes, "time" => time, "dual_obj" => dual_obj, 
                    "nbr_mtn" => nbr_mtn,  "nbr_sts_used" => nbr_sts_used, "status" => status, 
                    "mtn_stations_used" => mtn_stations_used)
        
        solution = Solution(obj_val, sx, sy, su, sv, sw, s_rho, s_lambda, s_phi, other_info)
        
        return solution
    else
        return Dict("status" => status)
    end
end
