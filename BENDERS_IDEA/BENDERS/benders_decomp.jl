

function solve_spdual(spdual_model, sx, sy, instance, A_theta2, A_theta3, eq_nodes)
    A_M, L_M = instance["A_M"], instance["L_M"]
    d, max_flt = instance["d"], instance["maximum_flying_time"] 
    a_nodes, f, V = instance["a_nodes"], instance["initial_flying_time"], instance["V"]
     
    b4, b5 = Dict{Tuple{String,String},Float64}(), Dict{Tuple{String,String},Float64}()
    b6, b7 = Dict{String,Float64}(), Dict{Tuple{String,String},Float64}()
    b8 = Dict{Tuple{String,String},Float64}()
    
    for a in Set(A_M)
        i,j = a
        b4[a] = max_flt * (sx[a]) - (max_flt - d[i]) * (1 - sy[j])
        b8[a] = d[j] - max_flt * (1 - sx[a]) - max_flt * sy[j]
    end
    for a in Set(A_theta2)
        i,j = a
        b5[a] = d[j] + (max_flt - d[i] - d[j]) * (1 - sx[a])
    end
    for j in Set(L_M)
        b6[j] = max_flt - (max_flt - d[j]) * sy[j]
    end
    for a in Set(A_theta3)
        i,j = a
        b7[a] = d[j] - max_flt * (1 - sx[a])
    end
    theta1, theta2, sigma, theta3 = spdual_model[:theta1], spdual_model[:theta2], spdual_model[:sigma], spdual_model[:theta3]
    theta4, pi1, alpha1, alpha2 = spdual_model[:theta4], spdual_model[:pi1], spdual_model[:alpha1], spdual_model[:alpha2]

    @objective(spdual_model, Max,   sum((max_flt*(sx[(i,j)]) - (max_flt - d[i])*(1 - sy[j])) * theta1[(i,j)] for (i,j) in A_M) + 
                                    sum((d[j] + (max_flt - d[i] - d[j])*(1 - sx[(i,j)])) * theta2[(i,j)] for (i,j) in A_theta2) +
                                    sum((max_flt - (max_flt - d[j])*sy[j]) * sigma[j] for j in L_M) +
                                    sum((d[j] - max_flt*(1 - sx[(i,j)])) * theta3[(i,j)] for (i,j) in A_theta3) +
                                    sum((d[j] - max_flt*(1 - sx[(i,j)]) - max_flt*sy[j]) * theta4[(i,j)] for (i,j) in A_M) + 
                                    sum(f[k] * pi1[k] for k in a_nodes if (k in eq_nodes)) +
                                    sum(d[j] * alpha1[j] for j in V) + 
                                    sum(max_flt * alpha2[j] for j in V))

    @time begin
        optimize!(spdual_model)
    end
    spdual_status = termination_status(spdual_model)
    if spdual_status == MOI.OPTIMAL
        return (
            spdual_model = spdual_model, status = spdual_status, obj = objective_value(spdual_model),
            val_t1 = JuMP.value.(theta1), val_t2 = JuMP.value.(theta2),
            val_s = JuMP.value.(sigma), val_t3 = JuMP.value.(theta3),
            val_t4 = JuMP.value.(theta4), val_p = JuMP.value.(pi1), 
            val_a1 = JuMP.value.(alpha1), val_a2 = JuMP.value.(alpha2),
            time = round(solve_time(spdual_model), digits = 6))
    elseif spdual_status == MOI.DUAL_INFEASIBLE
        return (
            spdual_model = spdual_model, status = spdual_status, obj = Inf, 
            ray_t1 = Dict(a => MOI.get(spdual_model, MOI.VariablePrimal(), theta1[a]) for a in A_M),
            ray_t2 = Dict(a => MOI.get(spdual_model, MOI.VariablePrimal(), theta2[a])  for a in A_theta2),
            ray_s = Dict(j => MOI.get(spdual_model, MOI.VariablePrimal(), sigma[j]) for j in L_M),
            ray_t3 = Dict(a => MOI.get(spdual_model, MOI.VariablePrimal(), theta3[a]) for a in A_theta3),
            ray_t4 = Dict(a => MOI.get(spdual_model, MOI.VariablePrimal(), theta4[a]) for a in A_M),
            pi1_ray = Dict(k => MOI.get(spdual_model, MOI.VariablePrimal(), pi1[k]) for k in eq_nodes),
            ray_a1 = Dict(j => MOI.get(spdual_model, MOI.VariablePrimal(), alpha1[j]) for j in V),
            ray_a2 = Dict(j =>  MOI.get(spdual_model, MOI.VariablePrimal(), alpha2[j]) for j in V), 
            time = round(solve_time(spdual_model), digits = 6))
    else
        return ()
    end
end 

function benders_decomp1(instance_file, nbr_thread, silent, time_limit)
    instance = build_graph(instance_file, false)
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
    f = instance["initial_flying_time"]
    h = instance["initial_takeoff"]
    g = instance["initial_flying_day"]
    b = instance["b"]
    DT, AT = instance["DT"], instance["AT"]
    a_nodes = instance["a_nodes"]
    fl_nodes = instance["fl_nodes"]
    pred_A, succ_A = Dict(), Dict()
    pred_A, succ_A = update_neighborhood(A, pred_A, succ_A)
    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(A_M, pred_A_M, succ_A_M)
    pred_A_M_bar, succ_A_M_bar = Dict(), Dict()
    pred_A_M_bar, succ_A_M_bar = update_neighborhood(A_M_bar, pred_A_M_bar, succ_A_M_bar)

    a_day = instance["a_day"]
    nbr_TP = instance["nbr_TP"]
    TP = 1:nbr_TP

    t_debut = Base.time()
    # ===================== Master problem =====================
    master_model = Model(optimizer_with_attributes(Gurobi.Optimizer, "Threads" => nbr_thread))
    set_optimizer_attribute(master_model, "OutputFlag", 1)
    set_optimizer_attribute(master_model, "TimeLimit", time_limit) 
    set_optimizer_attribute(master_model, "Presolve", 0)
    set_silent(master_model)


    # ===================== Decision variables =====================
    @variable(master_model, x[k in A], Bin)                          #Arc (i,j) selected
    @variable(master_model, y[j in L_M], Bin)                        #1 if a maintenance is performed at node j 
    @variable(master_model, Z >= 0)
    
    # ===================== Route Constraints =====================
    @constraint(master_model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(master_model, sum(x[i] for i in A_T) == nbr_K)
    @constraint(master_model, c1[i in fl_nodes], sum(x[(i,j)] for j in get(succ_A, i, [])) == 1)
    @constraint(master_model, c2[i in V_wt_st], sum(x[(j,i)] for j in get(pred_A, i, [])) == sum(x[(i,j)] for j in get(succ_A, i, [])))
    @constraint(master_model, c3[j in L_M], y[j] <= sum(x[(i, j)] for (i, j) in A_M))

    # ===================== Master Objective function =====================
    @objective(master_model, Min, Z)
    
    
    # =========================== Subproblem dual =====================
    A_theta2 = [(i,j) for (i,j) in A if j != "t"]
    A_theta3 = [(i,j) for (i,j) in A_M_bar if j != "t"]
    pred_A_theta2, succ_A_theta2 = Dict(), Dict()
    pred_A_theta2, succ_A_theta2 = update_neighborhood(A_theta2, pred_A_theta2, succ_A_theta2)
    pred_A_theta3, succ_A_theta3 = Dict(), Dict()
    pred_A_theta3, succ_A_theta3 = update_neighborhood(A_theta3, pred_A_theta3, succ_A_theta3)
    eq_nodes = union(["s"], a_nodes)

    spdual_model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(spdual_model, "OutputFlag", 1)  # passe à 0 pour silence
    set_optimizer_attribute(spdual_model, "InfUnbdInfo", 1)  # détecte infaisable vs non borné
    set_optimizer_attribute(spdual_model, "DualReductions", 0) # utile pour LP dual
    set_silent(spdual_model)
    # set_optimizer_attribute(spdual_model, "Threads", nbr_thread)  # si tu veux régler threads

    # --- variables duales (signes selon la dérivation) ---
    
    spdual_model[:theta1] = @variable(spdual_model, theta1[a in A_M] >= 0)             #Contrainte (c4)
    spdual_model[:theta2] = @variable(spdual_model, theta2[a in A_theta2] <= 0)        #Contrainte (c5) 
    spdual_model[:sigma] = @variable(spdual_model, sigma[j in L_M] <= 0)              #Contrainte (c6)
    spdual_model[:theta3] = @variable(spdual_model, theta3[a in A_theta3] >= 0)        #Contrainte (c7)
    spdual_model[:theta4] = @variable(spdual_model, theta4[a in A_M] >= 0)             #Contrainte (c8)
    spdual_model[:pi1] = @variable(spdual_model, pi1[k in eq_nodes])                #Contrainte u_k=f_k et u_s = 0 
    spdual_model[:alpha1] = @variable(spdual_model, alpha1[j in V] >= 0)               #Contrainte de borne inf u_j >= d_j
    spdual_model[:alpha2] = @variable(spdual_model, alpha2[j in V] <= 0)               #Contrainte de borne sup u_j <= max_flt

    # --- Contraintes duales ---
    # 1) pour chaque j in L_M (variable rho_j >= 0): sum_{i: (i,j) in A_M} theta1_{ij} <= 1
    @constraint(spdual_model, c1[i in L_M], sum(theta1[(k,i)] for k in get(pred_A_M, i, [])) <= 1)

    # 2) pour chaque k in V (variable u_k): combinaison des duales <= 0
    sum_theta1_out = Dict(k => sum(theta1[(k,j)] for j in get(succ_A_M, k, []); init=0.0) for k in V)
    sum_theta2_in  = Dict(k => sum(theta2[(i,k)] for i in get(pred_A_theta2, k, []); init=0.0) for k in V)
    sum_theta2_out = Dict(k => sum(theta2[(k,j)] for j in get(succ_A_theta2, k, []); init=0.0) for k in V)
    sum_theta3_in  = Dict(k => sum(theta3[(i,k)] for i in get(pred_A_theta3, k, []); init=0.0) for k in V)
    sum_theta3_out = Dict(k => sum(theta3[(k,j)] for j in get(succ_A_theta3, k, []); init=0.0) for k in V)
    sum_theta4_in  = Dict(k => sum(theta4[(i,k)] for i in get(pred_A_M, k, []); init=0.0) for k in V)
    sum_theta4_out = Dict(k => sum(theta4[(k,j)] for j in get(succ_A_M, k, []); init=0.0) for k in V)
    alpha_1_2_terms = Dict(k => alpha1[k] + alpha2[k] for k in V)
    sigma_term = Dict(k => ((k in L_M) ? sigma[k] : 0) for k in V)
    pi1_term    = Dict(k => ((k in eq_nodes) ? pi1[k] : 0) for k in V)


    @constraint(spdual_model, c2[k in setdiff(V, union(a_nodes, L_M, ["s"]))], sum_theta1_out[k] + sum_theta2_in[k] - sum_theta2_out[k] + 
                                                                        sum_theta3_in[k] - sum_theta3_out[k] + sum_theta4_in[k]  - 
                                                                        sum_theta4_out[k] + alpha_1_2_terms[k] <= 0)

    @constraint(spdual_model, c3[k in L_M],    sum_theta1_out[k] + sum_theta2_in[k] - sum_theta2_out[k] + sum_theta3_in[k] - 
                                        sum_theta3_out[k] + sum_theta4_in[k]  - sum_theta4_out[k] + alpha_1_2_terms[k] + 
                                        sigma_term[k] <= 0)

    @constraint(spdual_model, c4[k in union(a_nodes, ["s"])],  sum_theta1_out[k] + sum_theta2_in[k] - sum_theta2_out[k] + sum_theta3_in[k] - 
                                                        sum_theta3_out[k] + sum_theta4_in[k]  - sum_theta4_out[k] + alpha_1_2_terms[k] + 
                                                        pi1_term[k] <= 0)
 
    
    n = 0
    epsilon = 1E-4
    boucle = true
    optimal = false
    obj_val = 0
    while boucle
        println("\nIteration ", n)
        # Résolution du problème maître restreint
        @time begin
            optimize!(master_model)
        end
        master_time = round(solve_time(master_model), digits = 6)
        println("Temps du master : ", master_time)

        master_status = termination_status(master_model)
        if master_status == MOI.INFEASIBLE || master_status == MOI.INFEASIBLE_OR_UNBOUNDED
            println("Problème non réalisable")
            break
        end
        master_obj = objective_value(master_model)
        sx, sy, sZ = round.(JuMP.value.(x)), round.(JuMP.value.(y)), JuMP.value.(Z)
        println("Objectif du master :", master_obj)
        
        spd_result = solve_spdual(spdual_model, sx, sy, instance, A_theta2, A_theta3, eq_nodes)
        spdual_model = spd_result.spdual_model
        println("Temps du sp_dual = ", spd_result.time)
        if spd_result.status == MOI.INFEASIBLE_OR_UNBOUNDED
            println("Sous-problème dual non réalisable: problème non borné")
            break
        elseif spd_result.status == MOI.DUAL_INFEASIBLE
            println("Sous-problème infaisable: ajout d'une coupe de faisabilité")
            c13 = @constraint(master_model, sum(spd_result.ray_t1[(i,j)]*(max_flt * (sx[(i,j)]) - (max_flt - d[i]) * (1 - y[j])) for (i,j) in A_M) + 
                                sum(spd_result.ray_t2[(i,j)]*(d[j] + (max_flt - d[i] - d[j]) * (1 - x[(i,j)])) for (i,j) in A_theta2) + 
                                sum(spd_result.ray_s[j]*(max_flt - (max_flt - d[j]) * y[j]) for j in L_M) + 
                                sum(spd_result.ray_t3[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)])) for (i,j) in A_theta3) + 
                                sum(spd_result.ray_t4[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)]) - max_flt * y[j]) for (i,j) in A_M) + 
                                sum(spd_result.pi1_ray[k]*(f[k]) for k in a_nodes) + 
                                sum(spd_result.ray_a1[j]*(d[j]) for j in V) + 
                                sum(spd_result.ray_a2[j]*(max_flt) for j in V) <= 0)
            #println(c13)
        elseif spd_result.status == MOI.OPTIMAL
            print("Sous-problème dual a une solution bornée: ajout d'une coupe d'optimalité\n")
            println("Objectif du sp_dual :", spd_result.obj)
            c14 = @constraint(master_model,   sum(spd_result.val_t1[(i,j)]*(max_flt * (sx[(i,j)]) - (max_flt - d[i]) * (1 - y[j])) for (i,j) in A_M) + 
                                        sum(spd_result.val_t2[(i,j)]*(d[j] + (max_flt - d[i] - d[j]) * (1 - x[(i,j)])) for (i,j) in A_theta2) + 
                                        sum(spd_result.val_s[j]*(max_flt - (max_flt - d[j]) * y[j]) for j in L_M) + 
                                        sum(spd_result.val_t3[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)])) for (i,j) in A_theta3) + 
                                        sum(spd_result.val_t4[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)]) - max_flt * y[j]) for (i,j) in A_M) + 
                                        sum(spd_result.val_p[k]*(f[k]) for k in a_nodes) + 
                                        sum(spd_result.val_a1[j]*(d[j]) for j in V) + 
                                        sum(spd_result.val_a2[j]*(max_flt) for j in V) <= Z)  
            #println(c14)

        end
        if spd_result.obj - sZ <= epsilon
            boucle = false
            optimal = true
            obj_val = sZ
        end 
        n=n+1
    end
    t_ecoule = Base.time() - t_debut
    return Dict("obj" => obj_val, "time" => t_ecoule)
end 

