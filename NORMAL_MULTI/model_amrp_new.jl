using Dates

include("../structures.jl")

function model_amrp(env, instance_data, output_file, nbr_thread, silent, graph_reduc, time_limit, FH, FC, DY, MTN_CAP)
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

    #= if silent
        set_silent(model)
    end 
     =#
    #= gamma_f = 1       #1200 dollars par heure pour une maintenance
    gamma_t = 0      #2h(120 min) en moyenne par vol
    gamma_d = 0       #10h(600 min) en moyenne par jour =#
    
    ft_obj, tk_obj, fd_obj = 0, 0, 0

    # ===================== Decision variables =====================
    @variable(model, x[k in A], Bin)                          #Arc (i,j) selected
    @variable(model, y[j in L_M], Bin)                        #1 if a maintenance is performed at node j 

    # ===================== Route Constraints =====================
    @constraint(model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(model, sum(x[i] for i in A_T) == nbr_K)
    @constraint(model, c1[i in fl_nodes], sum(x[(i,j)] for j in get(successors, i, [])) == 1)
    @constraint(model, c2[i in V_wt_st], sum(x[(j,i)] for j in get(predecessors, i, [])) == sum(x[(i,j)] for j in get(successors, i, [])))
    @constraint(model, c3[j in L_M], y[j] <= sum(x[(i, j)] for i in get(pred_A_M, j, [])))
    @constraint(model, sum(y[j] for j in L_M) <= length(H_aircraft))
    
    
    #Valid inequalities
    #= @constraint(model, c32[j in V_wt_st], (length(get(pred_A_M_bar, j, []))-1) + sum(x[(i, j)] for i in get(pred_A_M_bar, j, [])) >= 1)
    @constraint(model, c33[i in V_wt_st], (length(get(succ_A_M_bar, i, []))-1) + sum(x[(i, j)] for j in get(succ_A_M_bar, i, [])) >= 1)
     =# 
    if FC
        # =========== Variables ===========                        
        @variable(model, tk[j] <= v[j in V] <= max_tk)                      #Accumulatad number of takeoff at node j
        @variable(model, lambda[j in L_M] >= 0, Int)                         #Remaining number of takeoff time at node j
    end 
    if FH 
        # =========== Variables ===========                        
        if graph_reduc 
            @variable(model, d_bar[j] <= u[j in V] <= F_bar[j])                  #Accumulatad flying time at node j 
        else
            @variable(model, d[j] <= u[j in V] <= max_flt)                  #Accumulatad flying time at node j 
        end 
        @variable(model, rho[j in L_M] >= 0)                              #Remaining flying time at node j

        # =========== Objective ===========
        gamma_f = 1     
        ft_obj = sum(gamma_f*rho[j] for j in L_M)

        # ========  Flying time constraints ========
        @constraint(model, c4[(i,j) in A_M], rho[j] - 76*lambda[j] >= max_flt*x[(i, j)] - u[i] - (max_flt - d[i])*(1 - y[j]))
        @constraint(model, [(i,j) in A; j != "t"], u[j] <= u[i] + d[j] + (max_flt - d[i] - d[j])*(1 - x[(i,j)]), base_name = "c5[($i,$j)]")
        @constraint(model, c6[j in L_M], u[j] <= max_flt -(max_flt - d[j])*y[j])
        @constraint(model, [(i,j) in A_M_bar; j != "t"], u[j] >= u[i] + d[j] - max_flt*(1 - x[(i,j)]), base_name = "c7[($i,$j)]")
        @constraint(model, c8[(i,j) in A_M], u[j] >= u[i] + d[j] - max_flt*(1 - x[(i,j)]) - max_flt*y[j])
        @constraint(model, u["s"] == 0)
        @constraint(model, c9[k in a_nodes], u[k] == f[k]) 
    end 

    if FC
        
        # =========== Objective ===========
        #= if tk_option == "MIN"
            gamma_t = other_data["fh_tk_min"]
        elseif tk_option == "MAX"
            gamma_t = other_data["fh_tk_max"]
        elseif tk_option == "MEDIAN"
            gamma_t = other_data["fh_tk_median"]
        elseif tk_option == "MEAN"
            gamma_t = other_data["fh_tk_mean"]
        else 
            gamma_t = 1
        end     =#
        #gamma_t = other_data["fh_tk"]
        gamma_t = 0
        tk_obj = sum(gamma_t*lambda[j] for j in L_M)

        # ========= Takeoff constraints ========
        @constraint(model, c10[(i,j) in A_M], lambda[j] >=  max_tk*x[(i, j)] - v[i] - (max_tk - tk[i])*(1 - y[j]))
        @constraint(model, [(i,j) in A; j != "t"], v[j] <= v[i] + tk[j] + (max_tk - tk[i] - tk[j])*(1 - x[(i,j)]), base_name = "c11[($i,$j)]")
        @constraint(model, c12[j in L_M], v[j] <= max_tk -(max_tk - tk[j])*y[j])
        @constraint(model, [(i,j) in A_M_bar; j != "t"], v[j] >= v[i] + tk[j] - max_tk*(1 - x[(i,j)]), base_name = "c13[($i,$j)]")
        @constraint(model, c14[(i,j) in A_M], v[j] >= v[i] + tk[j] - max_tk*(1 - x[(i,j)]) - max_tk*y[j])
        @constraint(model, v["s"] == 0)
        @constraint(model, c15[k in a_nodes], v[k] == h[k])    
    end
    if DY
        # =========== Variables ===========                        
        @variable(model, 0 <= w[j in V] <= max_day)                     #Accumulatad number of flying day at node j
        @variable(model, phi[j in L_M] >= 0, Int)                            #Remaining number of flying day at node j
        
        # =========== Objective ===========
        #= if fd_option == "MIN"
            gamma_d = instance["fh_day_min"]
        elseif fd_option == "MAX"
            gamma_d = instance["fh_day_max"]
        elseif fd_option == "MEDIAN"
            gamma_d = instance["fh_day_median"]
        elseif fd_option == "MEAN"
            gamma_d = instance["fh_day_mean"]
        else 
            gamma_d = 1
        end     =#
        
        gamma_d = 1

        #fd_obj = sum(gamma_d*phi[j] for j in L_M)

        # ========= Flying day constraints =========
        @constraint(model, c17[(i,j) in A_M], phi[j] >= max_day*x[(i, j)] - w[i] - (max_day - 1)*(1 - y[j]))
	    @constraint(model, [(i,j) in A; j != "t"], w[j] <= w[i] + b[(i,j)] + (max_day - b[(i,j)] - 1)*(1 - x[(i,j)]), base_name = "c18[($i,$j)]")
        @constraint(model, c19[j in L_M], w[j] <= max_day - (max_day - 1)*y[j])
        @constraint(model, [(i,j) in A_M_bar; j != "t"], w[j] >= w[i] + b[(i,j)] - max_day*(1 - x[(i,j)]), base_name = "c20[($i,$j)]")
        @constraint(model, c21[(i,j) in A_M], w[j] >= w[i] + b[(i,j)] - max_day*(1 - x[(i,j)]) - max_day*y[j])
        @constraint(model, w["s"] == 0)
        @constraint(model, c22[k in a_nodes], w[k] == g[k])
        @constraint(model, c23[j in V_wt_st], 1 <= w[j])      
    end

    if MTN_CAP
        @variable(model, z[m in MS, t in TP])                     #Number of maintenance opération at station m in period t
        #Maintenance capacity
        @constraint(model, c24[ms in MS, t in 1:nbr_TP], z[ms,t] == sum(y[i] for i in L_MS[ms] if fl_day[i] == t))
        @constraint(model, c25[ms in MS, t in 1:nbr_TP],  z[ms,t] <= ms_capacity[ms][t]) 
    end
    
    #@constraint(model, c28[j in L_M], fh_tk*(gamma_f-gamma_t-gamma_d)*lambda[j] <= (gamma_f-gamma_t-gamma_d)*rho[j])
    #= @constraint(model, c30[j in L_M], 3*(gamma_f+gamma_t-gamma_d)*phi[j] <= (gamma_f+gamma_t-gamma_d)*lambda[j])
    @constraint(model, c29[j in L_M], (gamma_f+gamma_t+gamma_d)*rho[j] <= 500*(gamma_f+gamma_t+gamma_d)*phi[j])
     =#

    # ===================== Objective general =====================
    @objective(model, Min, ft_obj + tk_obj + fd_obj)
    
    #write_to_file(model, "mon_modele.lp") 
    
    # ===================== Solve =====================

    optimize!(model)

    # ===================== Résultats =====================

    # Contrôle du statut
    status = termination_status(model)
    if status == MOI.OPTIMAL || status == MOI.FEASIBLE_POINT || status == MOI.TIME_LIMIT || status == MOI.INTERRUPTED
        sx, sy = Dict((i,j) => JuMP.value(x[(i,j)]) for (i,j) in A), Dict(j => JuMP.value(y[j]) for j in L_M)
        
        obj_val = objective_value(model)
        gap = round(relative_gap(model)*100, digits = 2)
        nbr_nodes =  MOI.get(model, MOI.NodeCount())
        dual_obj = round(objective_bound(model); digits = 2)
        
        time = round(solve_time(model), digits = 2)
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

        if FH
            su, s_rho = Dict(j => JuMP.value(u[j]) for j in V), Dict(j => round(Int, JuMP.value(rho[j])) for j in L_M)
            obj_rho = isempty(active_j) ? 0 : sum(s_rho[j] for j in active_j)
        else
            su, s_rho, obj_rho = Dict(j => 0 for j in V), Dict(j => 0 for j in L_M), 0
        end 
        if FC 
            sv, s_lambda = Dict(j => JuMP.value(v[j]) for j in V), Dict(j => round(Int, JuMP.value(lambda[j])) for j in L_M)
            obj_lambda = isempty(active_j) ? 0 : sum(s_lambda[j] for j in active_j)
        else
            sv, s_lambda, obj_lambda = Dict(j => 0 for j in V), Dict(j => 0 for j in L_M), 0
        end 
        if DY
            sw, s_phi = Dict(j => JuMP.value(w[j]) for j in V), Dict(j => round(Int, JuMP.value(phi[j])) for j in L_M)
            obj_phi = isempty(active_j) ? 0 : sum(s_phi[j] for j in active_j)
        else   
            sw, s_phi ,obj_phi = Dict(j => 0 for j in V), Dict(j => 0 for j in L_M), 0 
        end

        other_info = Dict("obj_rho" => obj_rho,"obj_lambda" => obj_lambda, "obj_phi" => obj_phi, 
                    "gap" => gap, "nbr_nodes" => nbr_nodes, "time" => time, "dual_obj" => dual_obj, 
                    "nbr_mtn" => nbr_mtn,  "nbr_sts_used" => nbr_sts_used, "status" => status, 
                    "mtn_stations_used" => mtn_stations_used)
        
        solution = Solution(obj_val, sx, sy, su, sv, sw, s_rho, s_lambda, s_phi, other_info)
        
        if MTN_CAP
            sz = JuMP.value.(z)
            output_results["z"] = sz
        end

        return solution
    else
        temp = Dict(j => 0 for j in V)
        other_info = Dict("status" => status)
        solution = 
        return Dict("status" => status)
    end
end
