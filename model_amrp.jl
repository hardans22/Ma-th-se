using Dates

include("build_graph.jl")

function model_amrp(instance_file, nbr_thread, silent, preprocess, time_limit)
    @time begin
        instance = build_graph(instance_file, preprocess)
    end
    #instance = preprocess(instance)
    nbr_K, nbr_FL = instance["number_of_aircrafts"], instance["number_of_flight_legs"]
    max_flt = instance["maximum_flying_time"]
    max_tk = instance["maximum_takeoff"]
    max_day = instance["maximum_flying_day"]
    mtn_stations, ms_capacity = instance["maintenance_stations"], instance["mtn_station_capacity"]
    d, tk = instance["d"], instance["tk"]
    A_S, A_K, A_F, A_T = instance["A_S"], instance["A_K"], instance["A_F"], instance["A_T"]
    A_M, A, L_M, A_M_bar = instance["A_M"], instance["A"], instance["L_M"], instance["A_M_bar"]
    V_wt_st, V, a_nodes = instance["V_wt_st"], instance["V"], instance["a_nodes"]
    MS, L_MS = instance["maintenance_stations"], instance["L_MS"]
    M_FL_O, M_FL_D = instance["M_FL_O"], instance["M_FL_D"]
    TRT = instance["turn_around_time"]
    f = instance["initial_flying_time"]
    h = instance["initial_takeoff"]
    g = instance["initial_flying_day"]
    b, b_bar = instance["b"], instance["b_bar"]
    b_bis = instance["b_bis"]
    a_nodes = instance["a_nodes"]
    fl_nodes = instance["fl_nodes"]
    predecessors, successors = Dict(), Dict()
    predecessors, successors = update_neighborhood(A, predecessors, successors)
    
    #= for i in V
        println("$i : ", b_bis[i])
    end
    =#
    a_day = instance["a_day"]
    nbr_TP = instance["nbr_TP"]
    TP = 1:nbr_TP
    #=exp_part = instance["exp_part"]
    init_level_ep = instance["init_level_ep"] 
    rate_ep = instance["rate_ep"]
    =#
    cs = 3 
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
    if silent
        set_silent(model)
    end 
    # ===================== Decision variables =====================
    @time begin
        @variable(model, x[k in A], Bin)                                #Arc (i,j) selected
        @variable(model, y[j in L_M], Bin)                              #Perform maintenance before j 
        @variable(model, d_bar[j] <= u[j in V] <= F_bar[j])                  #Accumulatad flying time at node j 
        #@variable(model, d[j] <= u[j in V] <= max_flt)                  #Accumulatad flying time at node j 
        @variable(model, tk[j] <= v[j in V] <= max_tk)                      #Accumulatad number of takeoff at node j
        @variable(model, b_bis[j] <= w[j in V] <= max_day)                     #Accumulatad number of flying day at node j
        @variable(model, rho[j in L_M] >= 0)                            #Remaining flying time at node j
        @variable(model, lambda[j in L_M] >= 0)                         #Remaining number of takeoff time at node j
        @variable(model, phi[j in L_M] >= 0)                            #Remaining number of flying day at node j
        
        @variable(model, z[m in MS, t in TP])                     #Number of maintenance opération at station m in period t
        #Définition provisoire
        MS_R = MS[1:2]
        MS_L = MS[3:end]
        SP = 1:3
        req = [rand(1:5) for _ in SP]
        D =Dict(ms => rand(0:5, length(SP), nbr_TP) for ms in MS)
        Q_max = rand(0:20, length(SP), nbr_TP)

        @variable(model, I[m in MS, k in SP, t in TP] >= 0)     #Inventory level of spare part p at station m in period t
        @variable(model, Q[m in MS, k in SP, t in TP] >= 0)     #Order quantity of spare part p at station m in period t
        @variable(model, S[m in MS_R, m1 in MS_L, k in SP, t in TP] >= 0)     #Transfert quantity of spare part p from station m to m1 in period t
        
        # ===================== Objective function =====================
        gamma_f = 20                #1200 dollars par heure pour une maintenance
        gamma_t = gamma_f*150       #2.5h(150 min) en moyenne par vol
        gamma_d = gamma_f*600       #10h(600 min) en moyenne par jour
        trt_wgt = gamma_f
        fl_wgt = 35*1000
        ft_obj = sum(gamma_f*rho[j] for j in L_M)
        tk_obj = sum(gamma_t*lambda[j] for j in L_M)
        fd_obj = sum(gamma_d*phi[j] for j in L_M)
        #trt_obj =  trt_wgt*sum((DT[j]-AT[i]-TRT)*x[(i,j)] for (i,j) in A_F) 
        
        @objective(model, Min, ft_obj + tk_obj + fd_obj)
        #@objective(model, Min, ft_obj + tk_obj + fd_obj + trt_obj)
        
        
        # ===================== Constraints =====================
        #@constraint(model, sum(x[i] for i in A_S) <= nbr_K)
        @constraint(model, sum(x[i] for i in A_S) == nbr_K)
        @constraint(model, sum(x[i] for i in A_T) == nbr_K)
        @constraint(model, c1[i in fl_nodes], sum(x[(i,j)] for j in get(successors, i, [])) == 1)
        @constraint(model, c2[i in V_wt_st], sum(x[(j,i)] for j in get(predecessors, i, [])) == sum(x[(i,j)] for j in get(successors, i, [])))
        @constraint(model, c3[j in L_M], y[j] <= sum(x[(i, j)] for (i, j) in A_M))
        
        #Flying time constraints
        @constraint(model, c4[(i,j) in A_M], rho[j] >= max_flt*x[(i, j)] - u[i] - (max_flt - d[i])*(1 - y[j]))
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
        #@constraint(model, c9[k in a_nodes], u[k] == f[k])
        
        #Takeoff constraints 
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
        #@constraint(model, c15[k in a_nodes], v[k] == h[k])
        #@constraint(model, c16[j in V_wt_st], tk[j] <= v[j])
        
        #Flying day constraints
        #= @constraint(model, c17[(i,j) in A_M], phi[j] >= max_day*x[(i, j)] - w[i] - (max_day - b_bis[i])*(1 - y[j]))
        for (i, j) in A
            if j != "t" && i != "s"
                @constraint(model, w[j] <= w[i] + b_bis[j]-1 + (max_day - b_bis[i] - b_bis[j])*(1 - x[(i,j)]), base_name = "c18[($i,$j)]")
            end
        end
        @constraint(model, c19[j in L_M], w[j] <= max_day - (max_day - b_bis[j])*y[j])
        for (i, j) in A_M_bar
            if j != "t" && i != "s"
                @constraint(model, w[j] >= w[i] + b_bis[j]-1 - max_day*(1 - x[(i,j)]), base_name = "c20[($i,$j)]")
            end
        end
        @constraint(model, c21[(i,j) in A_M], w[j] >= w[i] + b_bis[j]-1 - max_day*(1 - x[(i,j)]) - max_day*y[j])
        @constraint(model, w["s"] == 0)
        @constraint(model, c22[k in a_nodes], w[k] == g[k])
         =##@constraint(model, c23[j in V_wt_st], 1 <= w[j]) 
    
        #Maintenance capacity
        @constraint(model, c24[ms in MS, t in 1:nbr_TP], z[ms,t] == sum(y[i] for i in L_MS[ms] if a_day[i] == t))
        @constraint(model, c25[ms in MS, t in 1:nbr_TP],  z[ms,t] <= ms_capacity[ms][t])
        
        #Inventory constraints
        @constraint(model, c26[ms in MS_R, k in SP], I[ms,k,1] == 5)

        @constraint(model, c28[ms in MS_R, k in SP, t in 2:nbr_TP], 
                    I[ms,k,t] == I[ms,k,t-1] + Q[ms,k,t] - sum(S[ms,l,k,t] for l in MS_L) - req[k]*z[ms,t] - D[ms][k,t] )

        @constraint(model, c29[ms in MS_L, k in SP, t in 2:nbr_TP], 
                    I[ms,k,t] == I[ms,k,t-1] + Q[ms,k,t] + sum(S[l,ms,k,t] for l in MS_R) - req[k]*z[ms,t] - D[ms][k,t] )

        @constraint(model, c30[k in SP, t in TP], sum(Q[ms,k,t] for ms in MS) <= Q_max[k,t])

        #write_to_file(model, "model.lp")
        #= 
        @constraint(model, c25[ms in MS, t in 1:nbr_TP], z[ms,t] == sum(y[i]*a[(i,t)] for i in L_MS[ms]))

        @constraint(model, c26[p in exp_part, m in MS, t in 2:nbr_TP], se[p,m,t] == se[p,m,t-1] + qe[p,m,t] - rate_ep[p]*z[m,t] )
        @constraint(model, c27[p in exp_part, m in MS], se[p,m,1] == init_level_ep[p] + qe[p,m,1] - rate_ep[p]*z[m,1])
         =#
        
    end 
    optimize!(model)



    # ===================== Résultats =====================

    # Contrôle du statut
    status = termination_status(model)
    if status == MOI.OPTIMAL || status == MOI.FEASIBLE_POINT || status == MOI.TIME_LIMIT || status == MOI.INTERRUPTED
        obj_val = objective_value(model)
        sx, sy, su, sv, sw = JuMP.value.(x), JuMP.value.(y), JuMP.value.(u), JuMP.value.(v), JuMP.value.(w)
        sz, sI, sQ, sS = JuMP.value.(z), JuMP.value.(I), JuMP.value.(Q), JuMP.value.(S) 
        s_rho, s_lambda, s_phi = JuMP.value.(rho), JuMP.value.(lambda), JuMP.value.(phi)
        gap = round(relative_gap(model)*100, digits = 4)
        nbr_nodes =  MOI.get(model, MOI.NodeCount())
        dual_obj = objective_bound(model)
        time = round(solve_time(model), digits = 4)
        nbr_mtn = sum(sy)

        println("INVENTORY")
        for ms in MS
            println("Station : $ms")
            for k in SP
                for t in TP
                    if sI[ms,k,t] > 0
                        println("\nInventaire pour $k à $t: $ms\n")
                    end
                    if sQ[ms,k,t] > 0
                        println("\nQuantité commandée pour $k à $t: $ms\n")
                    end
                end
            end
        end

        println("\n\n\nQUANTITY TRANFERED ")
        for ms in MS_R
            for l in MS_L
                println("\nTransfert de $ms vers $l \n",sS[ms,l,:,:])
            end
        end




        mtn_stations_used = []
        active_j_set = Set(j for j in L_M if value(y[j]) == 1)
        remaining_stations = setdiff(Set(mtn_stations), mtn_stations_used)

        for ms in remaining_stations
            if !isempty(intersect(Set(L_MS[ms]), active_j_set))
                push!(mtn_stations_used, ms)
            end
        end
        
        nbr_sts_used = length(mtn_stations_used)

        return Dict("instance" => instance, "obj" => obj_val, "x" => sx, "y" => sy, "u" => su, "v" => sv, 
                    "w" => sw, "rho" => s_rho, "phi" => s_phi, "lambda" => s_lambda, "gap" => gap, 
                    "nbr_nodes" => nbr_nodes, "time" => time, "dual_obj" => dual_obj, "nbr_mtn" => nbr_mtn, 
                    "nbr_sts_used" => nbr_sts_used, "mtn_stations_used" => mtn_stations_used, 
                    "status" => status, "arc_reduc" => instance["arc_reduc"], "node_reduc" => instance["node_reduc"])
    else
        return Dict("status" => status)
    end
end