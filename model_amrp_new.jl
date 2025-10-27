using Dates

include("build_graph.jl")

function model_amrp(instance_file, nbr_thread, silent, preprocess, time_limit, FH, FC, DY, MTN_CAP)
    instance = build_graph(instance_file, preprocess)
    #instance = preprocess(instance)
    nbr_K = instance["number_of_aircrafts"]
    max_flt = instance["maximum_flying_time"]
    max_tk = instance["maximum_takeoff"]
    max_day = instance["maximum_flying_day"]
    mtn_stations, ms_capacity = instance["maintenance_stations"], instance["mtn_station_capacity"]
    d, tk = instance["d"], instance["tk"]
    A_S, A_T = instance["A_S"], instance["A_T"]
    A_M, A, L_M, A_M_bar = instance["A_M"], instance["A"], instance["L_M"], instance["A_M_bar"]
    V_wt_st, V, a_nodes = instance["V_wt_st"], instance["V"], instance["a_nodes"]
    MS, L_MS = instance["maintenance_stations"], instance["L_MS"]
    M_FL_O, M_FL_D = instance["M_FL_O"], instance["M_FL_D"]
    TRT = instance["turn_around_time"]
    f, h, g = instance["initial_flying_time"], instance["initial_takeoff"], instance["initial_flying_day"]
    b, DT, AT = instance["b"], instance["DT"], instance["AT"]
    a_nodes, fl_nodes= instance["a_nodes"], instance["fl_nodes"]
    fh_day, fh_tk = instance["fh_day"], instance["fh_tk"]
    
    predecessors, successors = Dict(), Dict()
    predecessors, successors = update_neighborhood(A, predecessors, successors)
    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(A, pred_A_M, succ_A_M)

    a_day = instance["a_day"]
    nbr_TP = instance["nbr_TP"]
    TP = 1:nbr_TP
    
    if preprocess
        F_bar = instance["F_bar"]
        d_bar = instance["d_bar"]
    else
        F_bar = Dict(j => instance["maximum_flying_time"] for j in V)
        d_bar = Dict(j => d[j] for j in V)
        instance["F_bar"], instance["d_bar"] = F_bar, d_bar
    end
    
    # ===================== Model =====================
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "Threads" => nbr_thread))
    set_optimizer_attribute(model, "OutputFlag", 1)
    set_optimizer_attribute(model, "TimeLimit", time_limit) 
    set_optimizer_attribute(model, "Presolve", 0)

    if silent
        set_silent(model)
    end 
    #= gamma_f = 0       #1200 dollars par heure pour une maintenance
    gamma_t = 1      #2h(120 min) en moyenne par vol
    gamma_d = 0       #10h(600 min) en moyenne par jour
    =#
    ft_obj, tk_obj, fd_obj = 0, 0, 0

    # ===================== Decision variables =====================
    @variable(model, x[k in A], Bin)                          #Arc (i,j) selected
    @variable(model, y[j in L_M], Bin)                        #1 if a maintenance is performed at node j 

    # ===================== Route Constraints =====================
    @constraint(model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(model, sum(x[i] for i in A_T) == nbr_K)
    @constraint(model, c1[i in fl_nodes], sum(x[(i,j)] for j in get(successors, i, [])) == 1)
    @constraint(model, c2[i in V_wt_st], sum(x[(j,i)] for j in get(predecessors, i, [])) == sum(x[(i,j)] for j in get(successors, i, [])))
    @constraint(model, c3[j in L_M], y[j] <= sum(x[(i, j)] for (i, j) in A_M))
    #@constraint(model, c35, sum(y[j] for j in L_M) <= 5)

    #Valid inequalities
    #= @constraint(model, c32[j in V_wt_st], (length(get(pred_A_M_bar, j, []))-1) + sum(x[(i, j)] for i in get(pred_A_M_bar, j, [])) >= 1)
    @constraint(model, c33[i in V_wt_st], (length(get(succ_A_M_bar, i, []))-1) + sum(x[(i, j)] for j in get(succ_A_M_bar, i, [])) >= 1)
     =# 
    #= for k in a_nodes 
        println()
        println()
        println()
        println(minimum(d[j] for j in get(succ_A_M, k, [])))
        println()
        println()
        println()
        if d[k] + minimum(d[j] for j in get(succ_A_M, k, [])) > max_flt
            @constraint(model, sum(y[j] for j in get(succ_A_M, k, [])) == 1)
        end
    end 
 =#
    if FH 
        # =========== Variables ===========                        
        @variable(model, d_bar[j] <= u[j in V] <= F_bar[j])                  #Accumulatad flying time at node j 
        #@variable(model, d[j] <= u[j in V] <= max_flt)                  #Accumulatad flying time at node j 
        @variable(model, rho[j in L_M] >= 0)                            #Remaining flying time at node j

        # =========== Objective ===========
        gamma_f = 1      
        ft_obj = sum(gamma_f*rho[j] for j in L_M)

        # ========  Flying time constraints ========
        @constraint(model, c4[(i,j) in A_M], rho[j] >= max_flt*x[(i, j)] - u[i] - (max_flt - d[i])*(1 - y[j]))
        #@constraint(model, c27, sum(rho[j] for j in L_M) <= Int(max_flt))

        for (i, j) in A
            if j != "t"
                @constraint(model, u[j] <= u[i] + d[j] + (max_flt - d[i] - d[j])*(1 - x[(i,j)]), base_name = "c5[($i,$j)]")
            end
        end
        @constraint(model, c6[j in L_M], u[j] <= max_flt -(max_flt - d[j])*y[j])
        for (i, j) in A_M_bar
            if j != "t"
                @constraint(model, u[j] >= u[i] + d[j] - max_flt*(1 - x[(i,j)]), base_name = "c7[($i,$j)]")
            end
        end
        
        @constraint(model, c8[(i,j) in A_M], u[j] >= u[i] + d[j] - max_flt*(1 - x[(i,j)]) - max_flt*y[j])
        @constraint(model, u["s"] == 0)
        @constraint(model, c9[k in a_nodes], u[k] == f[k])
        #@constraint(model, c35[j in L_M], (u[j]-d_bar[j])*(max_flt-u[j]) >= 0)  
          
    end 
    if FC
        # =========== Variables ===========                        
        @variable(model, tk[j] <= v[j in V] <= max_tk)                      #Accumulatad number of takeoff at node j
        @variable(model, lambda[j in L_M] >= 0, Int)                         #Remaining number of takeoff time at node j
        
        # =========== Objective ===========
        gamma_t = fh_tk     
        tk_obj = sum(gamma_t*lambda[j] for j in L_M)

        # ========= Takeoff constraints ========
        @constraint(model, c10[(i,j) in A_M], lambda[j] >= max_tk*x[(i, j)] - v[i] - (max_tk - tk[i])*(1 - y[j]))
            
        for (i, j) in A
            if j != "t" 
                @constraint(model, v[j] <= v[i] + tk[j] + (max_tk - tk[i] - tk[j])*(1 - x[(i,j)]), base_name = "c11[($i,$j)]")
            end
        end
        @constraint(model, c12[j in L_M], v[j] <= max_tk -(max_tk - tk[j])*y[j])
        for (i, j) in A_M_bar
            if j != "t" 
                @constraint(model, v[j] >= v[i] + tk[j] - max_tk*(1 - x[(i,j)]), base_name = "c13[($i,$j)]")
            end
        end
        @constraint(model, c14[(i,j) in A_M], v[j] >= v[i] + tk[j] - max_tk*(1 - x[(i,j)]) - max_tk*y[j])
        @constraint(model, v["s"] == 0)
        @constraint(model, c15[k in a_nodes], v[k] == h[k])    
    end
    if DY
        # =========== Variables ===========                        
        @variable(model, 0 <= w[j in V] <= max_day)                     #Accumulatad number of flying day at node j
        @variable(model, phi[j in L_M] >= 0, Int)                            #Remaining number of flying day at node j
        
        # =========== Objective ===========
        gamma_d = fh_day
        fd_obj = sum(gamma_d*phi[j] for j in L_M)

        # ========= Flying day constraints =========
        @constraint(model, c17[(i,j) in A_M], phi[j] >= max_day*x[(i, j)] - w[i] - (max_day - 1)*(1 - y[j]))
        #@constraint(model, c29, sum(phi[j] for j in L_M) <= Int(max_day))

	    for (i, j) in A
            if j != "t"
                @constraint(model, w[j] <= w[i] + b[(i,j)] + (max_day - b[(i,j)] - 1)*(1 - x[(i,j)]), base_name = "c18[($i,$j)]")
            end
        end 
        @constraint(model, c19[j in L_M], w[j] <= max_day - (max_day - 1)*y[j])
        for (i, j) in A_M_bar
            if j != "t" 
                @constraint(model, w[j] >= w[i] + b[(i,j)] - max_day*(1 - x[(i,j)]), base_name = "c20[($i,$j)]")
            end
        end
        @constraint(model, c21[(i,j) in A_M], w[j] >= w[i] + b[(i,j)] - max_day*(1 - x[(i,j)]) - max_day*y[j])
        @constraint(model, w["s"] == 0)
        @constraint(model, c22[k in a_nodes], w[k] == g[k])
        @constraint(model, c23[j in V_wt_st], 1 <= w[j])      
    end

    if MTN_CAP
        @variable(model, z[m in MS, t in TP])                     #Number of maintenance opération at station m in period t
        #Maintenance capacity
        @constraint(model, c24[ms in MS, t in 1:nbr_TP], z[ms,t] == sum(y[i] for i in L_MS[ms] if a_day[i] == t))
        @constraint(model, c25[ms in MS, t in 1:nbr_TP],  z[ms,t] <= ms_capacity[ms][t]) 
    end
    
    #= @constraint(model, c28[j in L_M], 175*(gamma_f-gamma_t-gamma_d)*lambda[j] <= (gamma_f-gamma_t-gamma_d)*rho[j])
    @constraint(model, c30[j in L_M], 3*(gamma_f+gamma_t-gamma_d)*phi[j] <= (gamma_f+gamma_t-gamma_d)*lambda[j])
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
        sx, sy = round.(JuMP.value.(x)), round.(Int,JuMP.value.(y))
        obj_val = objective_value(model)
        gap = round(relative_gap(model)*100, digits = 2)
        nbr_nodes =  MOI.get(model, MOI.NodeCount())
        dual_obj = round(objective_bound(model); digits = 2)
        time = round(solve_time(model), digits = 2)
        nbr_mtn = round(Int, sum(sy))

        mtn_stations_used = []
        active_j = Set(j for j in L_M if value(y[j]) >= 0.9)
        remaining_stations = setdiff(Set(mtn_stations), mtn_stations_used)

        for ms in remaining_stations
            if !isempty(intersect(Set(L_MS[ms]), active_j))
                push!(mtn_stations_used, ms)
            end
        end
        
        nbr_sts_used = round(Int, length(mtn_stations_used))

        if FH
            su, s_rho = JuMP.value.(u), round.(Int, JuMP.value.(rho))
            obj_rho = isempty(active_j) ? 0 : sum(s_rho[j] for j in active_j)
        else
            su, s_rho, obj_rho = Dict(j => 0 for j in V), Dict(j => 0 for j in L_M), 0
        end 
        if FC 
            sv, s_lambda = JuMP.value.(v), round.(Int,JuMP.value.(lambda))
            obj_lambda = isempty(active_j) ? 0 : sum(s_lambda[j] for j in active_j)
        else
            sv, s_lambda, obj_lambda = Dict(j => 0 for j in V), Dict(j => 0 for j in L_M), 0
        end 
        if DY
            sw, s_phi = JuMP.value.(w), round.(Int,JuMP.value.(phi))
            obj_phi = isempty(active_j) ? 0 : sum(s_phi[j] for j in active_j)
        else   
            sw, s_phi ,obj_phi = Dict(j => 0 for j in V), Dict(j => 0 for j in L_M), 0 
        end

        output_results = Dict("instance" => instance, "obj" => obj_val, "obj_rho" => obj_rho,"obj_lambda" => obj_lambda,
                    "obj_phi" => obj_phi, "x" => sx, "y" => sy, "u" => su, "v" => sv, "w" => sw, "rho" => s_rho, 
                    "phi" => s_phi, "lambda" => s_lambda, "gap" => gap, "nbr_nodes" => nbr_nodes, "time" => time, 
                    "dual_obj" => dual_obj, "nbr_mtn" => nbr_mtn,  "nbr_sts_used" => nbr_sts_used, 
                    "mtn_stations_used" => mtn_stations_used, "status" => status, 
                    "arc_reduc" => instance["arc_reduc"], "node_reduc" => instance["node_reduc"])
        
        if MTN_CAP
            sz = JuMP.value.(z)
            output_results["z"] = sz
        end

        return output_results
    else
        return Dict("status" => status)
    end
end
