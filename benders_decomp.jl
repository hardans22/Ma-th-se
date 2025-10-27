
function benders_decomp(instance_file, nbr_thread, silent, time_limit)
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

    # ===================== Solve master =====================
    master_status = termination_status(master_model)
    println("Master problem status: ", master_status)
    
    
    # =========================== Subproblem primal =====================
    sp_model = direct_model(Gurobi.Optimizer())
    set_attribute(sp_model, "InfUnbdInfo", 1)
    set_silent(sp_model)
    @variable(sp_model, 0 <= u[j in V])                  #Accualpha1latad flying time at node j 
    @variable(sp_model, rho[j in L_M] >= 0)                            #Remaining flying time at node j

    # =========== Objective ===========
    sigma_f = 1      
    ft_obj = sum(sigma_f*rho[j] for j in L_M)
    @objective(sp_model, Min, ft_obj)
    
    n = 0
    epsilon = 1E-4
    boucle = true
    optimal = false
    while boucle
        println("\nIteration ", n)
        # Résolution du problème maître restreint
        optimize!(master_model)
        master_status = termination_status(master_model)
        if master_status == MOI.INFEASIBLE || master_status == MOI.INFEASIBLE_OR_UNBOUNDED
            println("Problème non réalisable")
            break
        end
        master_obj = objective_value(master_model)
        sx, sy, sZ = round.(JuMP.value.(x)), round.(JuMP.value.(y)), JuMP.value.(Z)
        println("Objectif du master :", master_obj)
        
        # ========  Flying time constraints ========
        A_without_t = [a for a in A if a[2] != "t"]
        A_M_bar_without_t = [a for a in A_M_bar if a[2] != "t"]

        @constraint(sp_model, c4[(i,j) in A_M], rho[j] >= max_flt*sx[(i, j)] - u[i] - (max_flt - d[i])*(1 - sy[j]))
        @constraint(sp_model, c5[(i,j) in A_without_t], u[j] <= u[i] + d[j] + (max_flt - d[i] - d[j]) * (1 - sx[(i,j)]))
        @constraint(sp_model, c6[j in L_M], u[j] <= max_flt -(max_flt - d[j])*sy[j])
        @constraint(sp_model, c7[(i,j) in A_M_bar_without_t], u[j] >= u[i] + d[j] - max_flt*(1 - sx[(i,j)]))
        @constraint(sp_model, c8[(i,j) in A_M], u[j] >= u[i] + d[j] - max_flt*(1 - sx[(i,j)]) - max_flt*sy[j])
        @constraint(sp_model, c9[k in a_nodes], u[k] == f[k])
        @constraint(sp_model, c10, u["s"] == 0)
        @constraint(sp_model, c11[j in V], d[j] <= u[j])
        @constraint(sp_model, c12[j in V], u[j] <= max_flt)
        
        optimize!(sp_model)
        sp_obj = objective_value(sp_model)
        sp_status = termination_status(sp_model)
        if sp_status == MOI.INFEASIBLE_OR_UNBOUNDED
            println("Sous-problème dual non réalisable: problème non borné")
            break
        elseif sp_status == MOI.INFEASIBLE
            println("Sous-problème infaisable: ajout d'une coupe de faisabilité")
            R_C4 = Dict(a => MOI.get(sp_model, MOI.ConstraintDual(), c4[a]) for a in A_M)
            R_C5 = Dict(a => MOI.get(sp_model, MOI.ConstraintDual(), c5[a]) for a in A_without_t)
            R_C6 = Dict(j => MOI.get(sp_model, MOI.ConstraintDual(), c6[j]) for j in L_M)
            R_C7 = Dict(a => MOI.get(sp_model, MOI.ConstraintDual(), c7[j]) for a in A_M_bar_without_t)
            R_C8 = Dict(a => MOI.get(sp_model, MOI.ConstraintDual(), c8[a]) for a in A_M)
            R_C9 = Dict(k => MOI.get(sp_model, MOI.ConstraintDual(), c9[k]) for k in a_nodes)
            R_C11 = Dict(j => MOI.get(sp_model, MOI.ConstraintDual(), c11[j]) for j in V)
            R_C12 = Dict(j => MOI.get(sp_model, MOI.ConstraintDual(), c12[j]) for j in V)   
            c13 = @constraint(master_model, sum(R_C4[(i,j)]*(max_flt * (sx[(i,j)]) - (max_flt - d[i]) * (1 - y[j])) for (i,j) in A_M) + 
                                sum(R_C5[(i,j)]*(d[j] + (max_flt - d[i] - d[j]) * (1 - x[(i,j)])) for (i,j) in A_without_t) + 
                                sum(R_C6[j]*(max_flt - (max_flt - d[j]) * y[j]) for j in L_M) + 
                                sum(R_C7[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)])) for (i,j) in A_M_bar_without_t) + 
                                sum(R_C8[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)]) - max_flt * y[j]) for (i,j) in A_M) + 
                                sum(R_C9[k]*(f[k]) for k in a_nodes) + 
                                sum(R_C11[j]*(d[j]) for j in V) + 
                                sum(R_C12[j]*(max_flt) for j in V) <= 0)
            sp_obj = sZ+10
            #println(c13)
        else
            print("Sous-problème dual a une solution bornée: ajout d'une coupe d'optimalité")
            D_C4 = Dict(a => dual(c4[a]) for a in A_M)
            D_C5 = Dict(a => dual(c5[a]) for a in A_without_t)
            D_C6 = Dict(j => dual(c6[j]) for j in L_M)
            D_C7 = Dict(a => dual(c7[a]) for a in A_M_bar_without_t)
            D_C8 = Dict(a => dual(c8[a]) for a in A_M)
            D_C9 = Dict(k => dual(c9[k]) for k in a_nodes)
            D_C11 = Dict(j => dual(c11[j]) for j in V)
            D_C12 = Dict(j => dual(c12[j]) for j in V) 
            @constraint(master_model, sum(D_C4[(i,j)]*(max_flt * (sx[(i,j)]) - (max_flt - d[i]) * (1 - y[j])) for (i,j) in A_M) + 
                                sum(D_C5[(i,j)]*(d[j] + (max_flt - d[i] - d[j]) * (1 - x[(i,j)])) for (i,j) in A_without_t) + 
                                sum(D_C6[j]*(max_flt - (max_flt - d[j]) * y[j]) for j in L_M) + 
                                sum(D_C7[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)])) for (i,j) in A_M_bar_without_t) + 
                                sum(D_C8[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)]) - max_flt * y[j]) for (i,j) in A_M) + 
                                sum(D_C9[k]*(f[k]) for k in a_nodes) + 
                                sum(D_C11[j]*(d[j]) for j in V) + 
                                sum(D_C12[j]*(max_flt) for j in V) <= Z)  
        end
        if sp_obj - sZ <= epsilon
            boucle = false
            optimal = true
        end 
        n=n+1
        for constraint in c4
            delete(sp_model, constraint)
        end
        for constraint in c5
            delete(sp_model, constraint)
        end
        for constraint in c6
            delete(sp_model, constraint)
        end
        for constraint in c7
            delete(sp_model, constraint)
        end
        for constraint in c8
            delete(sp_model, constraint)
        end
        for constraint in c9
            delete(sp_model, constraint)
        end
        delete(sp_model, c10)  # c10 est une contrainte unique
        for constraint in c11
            delete(sp_model, constraint)
        end
        for constraint in c12
            delete(sp_model, constraint)
        end
        
        unregister(sp_model, :c4)
        unregister(sp_model, :c5)
        unregister(sp_model, :c6)
        unregister(sp_model, :c7)
        unregister(sp_model, :c8)
        unregister(sp_model, :c9)
        unregister(sp_model, :c10)
        unregister(sp_model, :c11)
        unregister(sp_model, :c12)
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
    
    @variable(spdual_model, theta1[a in A_M] >= 0)             #Contrainte (c4)
    @variable(spdual_model, theta2[a in A_theta2] <= 0)        #Contrainte (c5) 
    @variable(spdual_model, sigma[j in L_M] <= 0)              #Contrainte (c6)
    @variable(spdual_model, theta3[a in A_theta3] >= 0)        #Contrainte (c7)
    @variable(spdual_model, theta4[a in A_M] >= 0)             #Contrainte (c8)
    @variable(spdual_model, pi1[k in eq_nodes])                #Contrainte u_k=f_k et u_s = 0 
    @variable(spdual_model, alpha1[j in V] >= 0)               #Contrainte de borne inf u_j >= d_j
    @variable(spdual_model, alpha2[j in V] <= 0)               #Contrainte de borne sup u_j <= max_flt

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
        
        # ========  Flying time constraints ========
        b4, b5 = Dict{Tuple{String,String},Float64}(), Dict{Tuple{String,String},Float64}()
        b6, b7 = Dict{String,Float64}(), Dict{Tuple{String,String},Float64}()
        b8 = Dict{Tuple{String,String},Float64}()

        for a in A_M
            i,j = a
            b4[a] = max_flt * (sx[a]) - (max_flt - d[i]) * (1 - sy[j])
            b8[a] = d[j] - max_flt * (1 - sx[a]) - max_flt * sy[j]
        end
        for a in A_theta2
            i,j = a
            b5[a] = d[j] + (max_flt - d[i] - d[j]) * (1 - sx[a])
        end
        for j in L_M
            b6[j] = max_flt - (max_flt - d[j]) * sy[j]
        end
        for a in A_theta3
            i,j = a
            b7[a] = d[j] - max_flt * (1 - sx[a])
        end
        
        @objective(spdual_model, Max,   sum(b4[a] * theta1[a] for a in A_M) + sum(b5[a] * theta2[a] for a in A_theta2) +
                                        sum(b6[j] * sigma[j] for j in L_M) + sum(b7[a] * theta3[a] for a in A_theta3) +
                                        sum(b8[a] * theta4[a] for a in A_M) + sum(f[k] * pi1[k] for k in a_nodes if (k in eq_nodes)) +
                                        sum(d[j] * alpha1[j] for j in V) + sum(max_flt * alpha2[j] for j in V))

        @time begin
            optimize!(spdual_model)
        end
        sp_obj = objective_value(spdual_model)
        sp_status = termination_status(spdual_model)
        sp_time = round(solve_time(spdual_model), digits = 6)
        println("Temps du sp_dual : ", sp_time)
        
        if sp_status == MOI.INFEASIBLE_OR_UNBOUNDED
            println("Sous-problème dual non réalisable: problème non borné")
            break
        elseif sp_status == MOI.DUAL_INFEASIBLE
            println("Sous-problème infaisable: ajout d'une coupe de faisabilité")
            theta1_ray = Dict(a => MOI.get(spdual_model, MOI.VariablePrimal(), theta1[a]) for a in A_M)
            theta2_ray = Dict(a => MOI.get(spdual_model, MOI.VariablePrimal(), theta2[a])  for a in A_theta2)
            sigma_ray = Dict(j => MOI.get(spdual_model, MOI.VariablePrimal(), sigma[j]) for j in L_M)
            theta3_ray = Dict(a => MOI.get(spdual_model, MOI.VariablePrimal(), theta3[a]) for a in A_theta3)
            theta4_ray = Dict(a => MOI.get(spdual_model, MOI.VariablePrimal(), theta4[a]) for a in A_M)
            pi1_ray = Dict(k => MOI.get(spdual_model, MOI.VariablePrimal(), pi1[k]) for k in eq_nodes)
            alpha1_ray = Dict(j => MOI.get(spdual_model, MOI.VariablePrimal(), alpha1[j]) for j in V)
            alpha2_ray = Dict(j =>  MOI.get(spdual_model, MOI.VariablePrimal(), alpha2[j]) for j in V) 

            c13 = @constraint(master_model, sum(theta1_ray[(i,j)]*(max_flt * (sx[(i,j)]) - (max_flt - d[i]) * (1 - y[j])) for (i,j) in A_M) + 
                                sum(theta2_ray[(i,j)]*(d[j] + (max_flt - d[i] - d[j]) * (1 - x[(i,j)])) for (i,j) in A_theta2) + 
                                sum(sigma_ray[j]*(max_flt - (max_flt - d[j]) * y[j]) for j in L_M) + 
                                sum(theta3_ray[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)])) for (i,j) in A_theta3) + 
                                sum(theta4_ray[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)]) - max_flt * y[j]) for (i,j) in A_M) + 
                                sum(pi1_ray[k]*(f[k]) for k in a_nodes) + 
                                sum(alpha1_ray[j]*(d[j]) for j in V) + 
                                sum(alpha2_ray[j]*(max_flt) for j in V) <= 0)
            sp_obj = sZ+10
            #println(c13)
        else
            print("Sous-problème dual a une solution bornée: ajout d'une coupe d'optimalité\n")
            println("Objectif du sp_dual :", sp_obj)
            s_theta1, s_theta2, s_sigma, = JuMP.value.(theta1), JuMP.value.(theta2), JuMP.value.(sigma)
            s_theta3, s_theta4, s_pi1 = JuMP.value.(theta3), JuMP.value.(theta4), JuMP.value.(pi1)
            s_alpha1, s_alpha2 = JuMP.value.(alpha1), JuMP.value.(alpha2)
            c14 = @constraint(master_model,   sum(s_theta1[(i,j)]*(max_flt * (sx[(i,j)]) - (max_flt - d[i]) * (1 - y[j])) for (i,j) in A_M) + 
                                        sum(s_theta2[(i,j)]*(d[j] + (max_flt - d[i] - d[j]) * (1 - x[(i,j)])) for (i,j) in A_theta2) + 
                                        sum(s_sigma[j]*(max_flt - (max_flt - d[j]) * y[j]) for j in L_M) + 
                                        sum(s_theta3[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)])) for (i,j) in A_theta3) + 
                                        sum(s_theta4[(i,j)]*(d[j] - max_flt * (1 - x[(i,j)]) - max_flt * y[j]) for (i,j) in A_M) + 
                                        sum(s_pi1[k]*(f[k]) for k in a_nodes) + 
                                        sum(s_alpha1[j]*(d[j]) for j in V) + 
                                        sum(s_alpha2[j]*(max_flt) for j in V) <= Z)  
            #println(c14)

        end
        if sp_obj - sZ <= epsilon
            boucle = false
            optimal = true
            obj_val = sZ
        end 
        n=n+1
    end
    t_ecoule = Base.time() - t_debut
    return Dict("obj" => obj_val, "time" => t_ecoule)
end 

