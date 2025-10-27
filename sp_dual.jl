    #= # ----------------------------
    # Construire modèle dual
    # ----------------------------
    A_theta2 = [(i,j) for (i,j) in A if j != "t"]
    A_theta3 = [(i,j) for (i,j) in A_M_bar if j != "t"]
    pred_A_theta2, succ_A_theta2 = Dict(), Dict()
    pred_A_theta2, succ_A_theta2 = update_neighborhood(A_theta2, pred_A_theta2, succ_A_theta2)
    pred_A_theta3, succ_A_theta3 = Dict(), Dict()
    pred_A_theta3, succ_A_theta3 = update_neighborhood(A_theta3, pred_A_theta3, succ_A_theta3)
    
    eq_nodes = union(["s"], a_nodes)

    spdual_model = spdual_Model(Gurobi.Optimizer)
    set_optimizer_attribute(spdual_model, "OutputFlag", 1)  # passe à 0 pour silence
    set_optimizer_attribute(spdual_model, "InfUnbdInfo", 1)  # détecte infaisable vs non borné
    set_optimizer_attribute(spdual_model, "DualReductions", 0) # utile pour LP dual
    # set_optimizer_attribute(spdual_model, "Threads", nbr_thread)  # si tu veux régler threads

    # --- variables duales (signes selon la dérivation) ---
    
    @variable(spdual_model, theta1[a in A_M] >= 0)             #Contrainte (c4)
    @variable(spdual_model, theta2[a in A_theta2] <= 0)        #Contrainte (c5) 
    @variable(spdual_model, sigma[j in L_M] <= 0)              #Contrainte (c6)
    @variable(spdual_model, theta3[a in A_theta3] >= 0)        #Contrainte (c7)
    @variable(spdual_model, theta4[a in A_M] >= 0)             #Contrainte (c8)
    pi = Dict{String,JuMP.VariableRef}()
    for k in eq_nodes
        pi[k] = @variable(spdual_model, base_name = "pi_$k")   #Contrainte u_k=f_k et u_s = 0 
    end
    @variable(spdual_model, alpha1[j in V] >= 0)               #Contrainte de borne inf u_j >= d_j
    @variable(spdual_model, alpha2[j in V] <= 0)               #Contrainte de borne sup u_j <= max_flt

    # Calcul des termes de gauche pour les contraintes du probleme primal
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
    
    @objective(spdual_model, Max,  sum(b4[a] * theta1[a] for a in A_M) + sum(b5[a] * theta2[a] for a in A_theta2) +
                            sum(b6[j] * sigma[j] for j in L_M) + sum(b7[a] * theta3[a] for a in A_theta3) +
                            sum(b8[a] * theta4[a] for a in A_M) + sum(f[k] * pi[k] for k in a_nodes if haskey(pi,k)) +
                            sum(d[j] * alpha1[j] for j in V) + sum(max_flt * alpha2[j] for j in V))

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
    pi_term    = Dict(k => (haskey(pi, k) ? pi[k] : 0) for k in V)


    @constraint(spdual_model, c2[k in setdiff(V, union(a_nodes, L_M, ["s"]))], sum_theta1_out[k] + sum_theta2_in[k] - sum_theta2_out[k] + 
                                                                        sum_theta3_in[k] - sum_theta3_out[k] + sum_theta4_in[k]  - 
                                                                        sum_theta4_out[k] + alpha_1_2_terms[k] <= 0)

    @constraint(spdual_model,c_3[k in L_M],    sum_theta1_out[k] + sum_theta2_in[k] - sum_theta2_out[k] + sum_theta3_in[k] - 
                                        sum_theta3_out[k] + sum_theta4_in[k]  - sum_theta4_out[k] + alpha_1_2_terms[k] + 
                                        sigma_term[k] <= 0)

    @constraint(spdual_model, c4[k in union(a_nodes, ["s"])],  sum_theta1_out[k] + sum_theta2_in[k] - sum_theta2_out[k] + sum_theta3_in[k] - 
                                                        sum_theta3_out[k] + sum_theta4_in[k]  - sum_theta4_out[k] + alpha_1_2_terms[k] + 
                                                        pi_term[k] <= 0)
 
    # ----------------------------
    # Résolution
    # ----------------------------
    optimize!(spdual_model)

    status = termination_status(spdual_model)
    #obj = objective_value(spdual_model)
    println("Statut de terminaison : ", status)
    println("dual de terminaison : ", primal_status(spdual_model))
    @show shadow_price(c4["s"])
    =#
    #println("Valeur objective duale : ", obj)

    # --- Afficher quelques variables duales pour vérification ---
    #= println("\nQuelques variables duales (extrait) :")
    for a in A_M
        println("theta1[", a, "] = ", value(theta1[a]))
    end
    for a in A_theta2
        println("theta2[", a, "] = ", value(theta2[a]))
    end
    for j in L_M
        println("sigma[", j, "] = ", value(sigma[j]))
    end
    for j in V
        println("alpha1[", j, "] = ", value(alpha1[j]), "  alpha2[", j, "] = ", value(alpha2[j]))
    end
    for k in eq_nodes
        println("pi[", k, "] = ", value(pi[k]))
    end

 =#
