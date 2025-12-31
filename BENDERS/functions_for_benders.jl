function build_sp(env, instance)
    d, V = instance["d"], instance["V"]
    A, A_M, A_M_bar = instance["A"], instance["A_M"], instance["A_M_bar"]
    L_M = instance["L_M"]
    max_flt = instance["maximum_flying_time"]
    a_nodes = instance["a_nodes"]
    f = instance["initial_flying_time"]

    sp_model = Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(sp_model, "OutputFlag", 0)
    set_optimizer_attribute(sp_model, "InfUnbdInfo", 1)
    set_optimizer_attribute(sp_model, "DualReductions", 0)
    @variable(sp_model, x_copy[k in A], Bin)                          #Arc (i,j) selected
    @variable(sp_model, y_copy[j in L_M], Bin)        
    @variable(sp_model, d[j] <= u[j in V] <= max_flt)                #Accumulatad flying time at node j 
    #@variable(model, d[j] <= u[j in V] <= max_flt)                  #Accumulatad flying time at node j 
    @variable(sp_model, rho[j in L_M] >= 0)                            #Remaining flying time at node j

    # =========== Objective ===========
    @objective(sp_model, Min,  sum(rho[j] for j in L_M))

    # ========  Flying time constraints ========
    @constraint(sp_model, c4[(i,j) in A_M], rho[j] >= max_flt*x_copy[(i, j)] - u[i] - (max_flt - d[i])*(1 - y_copy[j]))

    for (i, j) in A
        if j != "t"
            @constraint(sp_model, u[j] <= u[i] + d[j] + (max_flt - d[i] - d[j])*(1 - x_copy[(i,j)]), base_name = "c5[($i,$j)]")
        end
    end
    @constraint(sp_model, c6[j in L_M], u[j] <= max_flt -(max_flt - d[j])*y_copy[j])
    for (i, j) in A_M_bar
        if j != "t"
            @constraint(sp_model, u[j] >= u[i] + d[j] - max_flt*(1 - x_copy[(i,j)]), base_name = "c7[($i,$j)]")
        end
    end
    
    @constraint(sp_model, c8[(i,j) in A_M], u[j] >= u[i] + d[j] - max_flt*(1 - x_copy[(i,j)]) - max_flt*y_copy[j])
    @constraint(sp_model, u["s"] == 0)
    @constraint(sp_model, c9[k in a_nodes], u[k] == f[k])

    return sp_model
end

function solve_sp(sp_model, x_val, y_val)
    # Fixer x_copy
    for k in keys(x_val)
        fix(sp_model[:x_copy][k], x_val[k])
    end

    # Fixer y_copy
    for j in keys(y_val)
        fix(sp_model[:y_copy][j], y_val[j])
    end
    optimize!(sp_model)
    sp_status = termination_status(sp_model)
    if sp_status == MOI.OPTIMAL
        return (status = sp_status, obj = objective_value(sp_model), time = round(solve_time(sp_model), digits = 4))
    else 
        #= if termination_status(sp_model) == MOI.INFEASIBLE
            # Calculer le conflit (IIS)
            compute_conflict!(sp_model)
            
            # Afficher les contraintes en conflit
            for (F, S) in list_of_constraint_types(sp_model)
                for con in all_constraints(sp_model, F, S)
                    if MOI.get(sp_model, MOI.ConstraintConflictStatus(), con) == MOI.IN_CONFLICT
                        println("Contrainte en conflit : ", con)
                        println("  ", constraint_object(con))
                    end
                end
            end
        end =#
        return (status = sp_status, obj = Inf)
    end
end 


function build_spdual(env, instance)
    A_M, L_M = instance["A_M"], instance["L_M"]
    #d, max_flt = instance["d"], instance["maximum_flying_time"] 
    eq_nodes, f, V = instance["eq_nodes"], instance["initial_flying_time"], instance["V"]
    # ===================== Subproblem dual =====================
    A_theta2 = instance["A_theta2"]
    A_theta3 = instance["A_theta3"]
    pred_A_theta2, succ_A_theta2 = Dict(), Dict()
    pred_A_theta2, succ_A_theta2 = update_neighborhood(A_theta2, pred_A_theta2, succ_A_theta2)
    pred_A_theta3, succ_A_theta3 = Dict(), Dict()
    pred_A_theta3, succ_A_theta3 = update_neighborhood(A_theta3, pred_A_theta3, succ_A_theta3)
    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(A_M, pred_A_M, succ_A_M)
    

    spdual_model = Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(spdual_model, "OutputFlag", 0)
    set_optimizer_attribute(spdual_model, "InfUnbdInfo", 1)
    set_optimizer_attribute(spdual_model, "DualReductions", 0)
    set_optimizer_attribute(spdual_model, "Threads", nbr_thread)
    # Variables duales
    @variable(spdual_model, theta1[a in A_M] >= 0)
    @variable(spdual_model, theta2[a in A_theta2] <= 0)
    @variable(spdual_model, sigma[j in L_M] <= 0)
    @variable(spdual_model, theta3[a in A_theta3] >= 0)
    @variable(spdual_model, theta4[a in A_M] >= 0)
    @variable(spdual_model, pi1[k in eq_nodes])
    @variable(spdual_model, alpha1[j in V] >= 0)
    @variable(spdual_model, alpha2[j in V] <= 0)

    # Contraintes duales
    @constraint(spdual_model, spd_c1[k in L_M], sum(theta1[(i,k)] for i in get(pred_A_M, k, [])) <= 1)

    sum_theta1_out = Dict(k => sum(theta1[(k,j)] for j in get(succ_A_M, k, []); init=0.0) for k in V)
    sum_theta2_in  = Dict(k => sum(theta2[(i,k)] for i in get(pred_A_theta2, k, []); init=0.0) for k in V)
    sum_theta2_out = Dict(k => sum(theta2[(k,j)] for j in get(succ_A_theta2, k, []); init=0.0) for k in V)
    sum_theta3_in  = Dict(k => sum(theta3[(i,k)] for i in get(pred_A_theta3, k, []); init=0.0) for k in V)
    sum_theta3_out = Dict(k => sum(theta3[(k,j)] for j in get(succ_A_theta3, k, []); init=0.0) for k in V)
    sum_theta4_in  = Dict(k => sum(theta4[(i,k)] for i in get(pred_A_M, k, []); init=0.0) for k in V)
    sum_theta4_out = Dict(k => sum(theta4[(k,j)] for j in get(succ_A_M, k, []); init=0.0) for k in V)
    

    @constraint(spdual_model, spd_c2[k in setdiff(V, union(eq_nodes, L_M))], 
        sum_theta1_out[k] + sum_theta2_in[k] - sum_theta2_out[k] + 
        sum_theta3_in[k] - sum_theta3_out[k] + sum_theta4_in[k] - 
        sum_theta4_out[k] + alpha1[k] + alpha2[k] <= 0)

    @constraint(spdual_model, spd_c3[k in L_M], 
        sum_theta1_out[k] + sum_theta2_in[k] - sum_theta2_out[k] + 
        sum_theta3_in[k] - sum_theta3_out[k] + sum_theta4_in[k] - 
        sum_theta4_out[k] + alpha1[k] + alpha2[k] + sigma[k] <= 0)

    @constraint(spdual_model, spd_c4[k in eq_nodes], 
        sum_theta1_out[k] + sum_theta2_in[k] - sum_theta2_out[k] + 
        sum_theta3_in[k] - sum_theta3_out[k] + sum_theta4_in[k] - 
        sum_theta4_out[k] + alpha1[k] + alpha2[k] + pi1[k] <= 0)

    # Enregistrer les variables dans le modèle
    spdual_model[:theta1], spdual_model[:theta2], spdual_model[:sigma], spdual_model[:theta3] = theta1, theta2, sigma, theta3
    spdual_model[:theta4], spdual_model[:pi1], spdual_model[:alpha1], spdual_model[:alpha2]  = theta4, pi1, alpha1, alpha2
    spdual_model[:spd_c1], spdual_model[:spd_c2], spdual_model[:spd_c3], spdual_model[:spd_c4] = spd_c1, spd_c2, spd_c3, spd_c4
    
    return spdual_model

end 

function solve_spdual(x_val, y_val, instance, spdual_model)
    A_M, L_M, A_theta2, A_theta3 = instance["A_M"], instance["L_M"], instance["A_theta2"], instance["A_theta3"]
    d, max_flt, = instance["d"], instance["maximum_flying_time"]
    a_nodes, f, V = instance["a_nodes"], instance["initial_flying_time"], instance["V"]
    
    # Récupérer les variables du modèle dual
    theta1, theta2, sigma = spdual_model[:theta1], spdual_model[:theta2], spdual_model[:sigma]
    theta3, theta4, pi1 = spdual_model[:theta3], spdual_model[:theta4], spdual_model[:pi1]
    alpha1, alpha2 = spdual_model[:alpha1], spdual_model[:alpha2]
    
    # Mettre à jour l'objectif avec les valeurs courantes de x et y
    @objective(spdual_model, Max,   
        sum((max_flt*x_val[(i,j)] - (max_flt - d[i])*(1 - sum(y_val[j,k] for k in a_nodes))) * theta1[(i,j)] for (i,j) in A_M) + 
        sum((d[j] + (max_flt - d[i] - d[j])*(1 - x_val[(i,j)])) * theta2[(i,j)] for (i,j) in A_theta2) +
        sum((max_flt - (max_flt - d[j])*sum(y_val[j,k] for k in a_nodes)) * sigma[j] for j in L_M) +
        sum((d[j] - max_flt*(1 - x_val[(i,j)])) * theta3[(i,j)] for (i,j) in A_theta3) +
        sum((d[j] - max_flt*(1 - x_val[(i,j)]) - max_flt*sum(y_val[j,k] for k in a_nodes)) * theta4[(i,j)] for (i,j) in A_M) + 
        sum(f[k] * pi1[k] for k in a_nodes) +
        sum(d[j] * alpha1[j] for j in V) + 
        sum(max_flt * alpha2[j] for j in V))

    optimize!(spdual_model)
    
    spdual_status = termination_status(spdual_model)
    
    if spdual_status == MOI.OPTIMAL
        #= for v in all_variables(spdual_model)
            rc = reduced_cost(v)
            if rc == 0
                println(rc)
                println("$(v) peut varier sans changer l'objectif")
            end
        end =#
        return (status = spdual_status, obj = objective_value(spdual_model), 
                val_t1 = value.(theta1), val_t2 = value.(theta2), val_s = value.(sigma), 
                val_t3 = value.(theta3), val_t4 = value.(theta4), val_p = value.(pi1), 
                val_a1 = value.(alpha1), val_a2 = value.(alpha2), 
                time = round(solve_time(spdual_model), digits = 4))
    elseif spdual_status == MOI.DUAL_INFEASIBLE
        return (status = spdual_status, obj = Inf, time = round(solve_time(spdual_model), digits = 4),
                ray_t1 = value.(theta1), ray_t2 = value.(theta2), ray_s = value.(sigma), 
                ray_t3 = value.(theta3), ray_t4 = value.(theta4), ray_p = value.(pi1), 
                ray_a1 = value.(alpha1), ray_a2 = value.(alpha2))
    else
        return (status = spdual_status, obj = nothing)
    end
end 

function solve_aux_spdual(x_val, y_val, core_x_val, core_y_val, instance, aux_spdual, sp_obj)
    A_M, L_M, A_theta2, A_theta3 = instance["A_M"], instance["L_M"], instance["A_theta2"], instance["A_theta3"]
    d, max_flt, eq_nodes = instance["d"], instance["maximum_flying_time"], instance["eq_nodes"] 
    a_nodes, f, V = instance["a_nodes"], instance["initial_flying_time"], instance["V"]
    
    # Récupérer les variables du modèle dual
    theta1, theta2, sigma = aux_spdual[:theta1], aux_spdual[:theta2], aux_spdual[:sigma]
    theta3, theta4, pi1 = aux_spdual[:theta3], aux_spdual[:theta4], aux_spdual[:pi1]
    alpha1, alpha2 = aux_spdual[:alpha1], aux_spdual[:alpha2]
    
    # Mettre à jour l'objectif avec les valeurs courantes de x et y
    @objective(aux_spdual, Max,   
        sum((max_flt*core_x_val[(i,j)] - (max_flt - d[i])*(1 - sum(core_y_val[j,k] for k in a_nodes))) * theta1[(i,j)] for (i,j) in A_M) + 
        sum((d[j] + (max_flt - d[i] - d[j])*(1 - core_x_val[(i,j)])) * theta2[(i,j)] for (i,j) in A_theta2) +
        sum((max_flt - (max_flt - d[j])*sum(core_y_val[j,k] for k in a_nodes)) * sigma[j] for j in L_M) +
        sum((d[j] - max_flt*(1 - core_x_val[(i,j)])) * theta3[(i,j)] for (i,j) in A_theta3) +
        sum((d[j] - max_flt*(1 - core_x_val[(i,j)]) - max_flt*sum(core_y_val[j,k] for k in a_nodes)) * theta4[(i,j)] for (i,j) in A_M) + 
        sum(f[k] * pi1[k] for k in a_nodes) +
        sum(d[j] * alpha1[j] for j in V) + 
        sum(max_flt * alpha2[j] for j in V))

    if aux_spdual.ext[:c5_ref] !== nothing
        delete(aux_spdual, aux_spdual.ext[:c5_ref])
    end

    aux_spdual.ext[:c5_ref] = @constraint(aux_spdual,
        sum((max_flt*x_val[(i,j)] - (max_flt - d[i])*(1 - sum(y_val[j,k] for k in a_nodes))) * theta1[(i,j)] for (i,j) in A_M) + 
        sum((d[j] + (max_flt - d[i] - d[j])*(1 - x_val[(i,j)])) * theta2[(i,j)] for (i,j) in A_theta2) +
        sum((max_flt - (max_flt - d[j])*sum(y_val[j,k] for k in a_nodes)) * sigma[j] for j in L_M) +
        sum((d[j] - max_flt*(1 - x_val[(i,j)])) * theta3[(i,j)] for (i,j) in A_theta3) +
        sum((d[j] - max_flt*(1 - x_val[(i,j)]) - max_flt*sum(y_val[j,k] for k in a_nodes)) * theta4[(i,j)] for (i,j) in A_M) + 
        sum(f[k] * pi1[k] for k in a_nodes if (k in eq_nodes)) +
        sum(d[j] * alpha1[j] for j in V) + 
        sum(max_flt * alpha2[j] for j in V) == sp_obj)

    optimize!(aux_spdual)
    spdual_status = termination_status(aux_spdual)
    
    if spdual_status == MOI.OPTIMAL
        return (status = spdual_status, obj = objective_value(aux_spdual), 
                val_t1 = value.(theta1), val_t2 = value.(theta2), val_s = value.(sigma), 
                val_t3 = value.(theta3), val_t4 = value.(theta4), val_p = value.(pi1), 
                val_a1 = value.(alpha1), val_a2 = value.(alpha2), 
                time = round(solve_time(aux_spdual), digits = 6))
    else
        return (status = spdual_status, obj = nothing)
    end
end 

function build_opti_cut(result, x, y, Z, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
    coupe = @build_constraint( sum(result.val_t1[(i,j)] * (max_flt * x[(i,j)] - (max_flt - d[i]) * (1 - sum(y[j,k] for k in a_nodes))) for (i,j) in A_M) + 
            sum(result.val_t2[(i,j)] * (d[j] + (max_flt - d[i] - d[j]) * (1 - x[(i,j)])) for (i,j) in A_theta2) + 
            sum(result.val_s[j] * (max_flt - (max_flt - d[j]) * sum(y[j,k] for k in a_nodes)) for j in L_M) + 
            sum(result.val_t3[(i,j)] * (d[j] - max_flt * (1 - x[(i,j)])) for (i,j) in A_theta3) + 
            sum(result.val_t4[(i,j)] * (d[j] - max_flt * (1 - x[(i,j)]) - max_flt * sum(y[j,k] for k in a_nodes)) for (i,j) in A_M) + 
            sum(result.val_p[k] * f[k] for k in a_nodes) + 
            sum(result.val_a1[j] * d[j] for j in V) + 
            sum(result.val_a2[j] * max_flt for j in V) <= Z)
    return coupe              
end

function build_fais_cut(result, x, y, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
    coupe = @build_constraint( sum(result.ray_t1[(i,j)] * (max_flt * x[(i,j)] - (max_flt - d[i]) * (1 - sum(y[j,k] for k in a_nodes))) for (i,j) in A_M) + 
        sum(result.ray_t2[(i,j)] * (d[j] + (max_flt - d[i] - d[j]) * (1 - x[(i,j)])) for (i,j) in A_theta2) + 
        sum(result.ray_s[j] * (max_flt - (max_flt - d[j]) * sum(y[j,k] for k in a_nodes)) for j in L_M) + 
        sum(result.ray_t3[(i,j)] * (d[j] - max_flt * (1 - x[(i,j)])) for (i,j) in A_theta3) + 
        sum(result.ray_t4[(i,j)] * (d[j] - max_flt * (1 - x[(i,j)]) - max_flt * sum(y[j,k] for k in a_nodes)) for (i,j) in A_M) + 
        sum(result.ray_p[k] * f[k] for k in a_nodes) + 
        sum(result.ray_a1[j] * d[j] for j in V) + 
        sum(result.ray_a2[j] * max_flt for j in V) <= 0)
    
end 

function compare_val(result, core_x_val, core_y_val, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
    
    val = sum(result.val_t1[(i,j)] * (max_flt * core_x_val[(i,j)] - (max_flt - d[i]) * (1 - core_y_val[j])) for (i,j) in A_M) + 
                    sum(result.val_t2[(i,j)] * (d[j] + (max_flt - d[i] - d[j]) * (1 - core_x_val[(i,j)])) for (i,j) in A_theta2) + 
                    sum(result.val_s[j] * (max_flt - (max_flt - d[j]) * core_y_val[j]) for j in L_M) + 
                    sum(result.val_t3[(i,j)] * (d[j] - max_flt * (1 - core_x_val[(i,j)])) for (i,j) in A_theta3) + 
                    sum(result.val_t4[(i,j)] * (d[j] - max_flt * (1 - core_x_val[(i,j)]) - max_flt * core_y_val[j]) for (i,j) in A_M) + 
                    sum(result.val_p[k] * f[k] for k in a_nodes) + 
                    sum(result.val_a1[j] * d[j] for j in V) + 
                    sum(result.val_a2[j] * max_flt for j in V)
    return val
end

function cut_comparison(result, result_aux, x_val, y_val, instance, A_theta2, A_theta3)
    A_M, L_M = instance["A_M"], instance["L_M"]
    d, max_flt = instance["d"], instance["maximum_flying_time"] 
    a_nodes, f, V = instance["a_nodes"], instance["initial_flying_time"], instance["V"]
    
    val = sum(result.val_t1[(i,j)] * (max_flt * x_val[(i,j)] - (max_flt - d[i]) * (1 - y_val[j])) for (i,j) in A_M) + 
                    sum(result.val_t2[(i,j)] * (d[j] + (max_flt - d[i] - d[j]) * (1 - x_val[(i,j)])) for (i,j) in A_theta2) + 
                    sum(result.val_s[j] * (max_flt - (max_flt - d[j]) * y_val[j]) for j in L_M) + 
                    sum(result.val_t3[(i,j)] * (d[j] - max_flt * (1 - x_val[(i,j)])) for (i,j) in A_theta3) + 
                    sum(result.val_t4[(i,j)] * (d[j] - max_flt * (1 - x_val[(i,j)]) - max_flt * y_val[j]) for (i,j) in A_M) + 
                    sum(result.val_p[k] * f[k] for k in a_nodes) + 
                    sum(result.val_a1[j] * d[j] for j in V) + 
                    sum(result.val_a2[j] * max_flt for j in V)

    val_aux = sum(result_aux.val_t1[(i,j)] * (max_flt * x_val[(i,j)] - (max_flt - d[i]) * (1 - y_val[j])) for (i,j) in A_M) + 
                    sum(result_aux.val_t2[(i,j)] * (d[j] + (max_flt - d[i] - d[j]) * (1 - x_val[(i,j)])) for (i,j) in A_theta2) + 
                    sum(result_aux.val_s[j] * (max_flt - (max_flt - d[j]) * y_val[j]) for j in L_M) + 
                    sum(result_aux.val_t3[(i,j)] * (d[j] - max_flt * (1 - x_val[(i,j)])) for (i,j) in A_theta3) + 
                    sum(result_aux.val_t4[(i,j)] * (d[j] - max_flt * (1 - x_val[(i,j)]) - max_flt * y_val[j]) for (i,j) in A_M) + 
                    sum(result_aux.val_p[k] * f[k] for k in a_nodes) + 
                    sum(result_aux.val_a1[j] * d[j] for j in V) + 
                    sum(result_aux.val_a2[j] * max_flt for j in V)
    println("VAL CUT NORMAL: $val")
    println("VAL CUT PARETO: $val_aux")
    #= if val >= val_aux
        return result
    else 
        #println("PARETOOOOOOOOO")
        return result_aux
    end =# 
end 

function y_neighborhood(x_val, y_val, acritik_set, A, L_M, A_M)
    y_list = []
    aircraft_paths = build_all_paths(x_val, A)
    
    for ac in acritik_set
        path = aircraft_paths[ac]
        len_p = length(path)
        mtn_node = []
        for i in 2:len_p  # Commencer à 2 car on a besoin d'un prédécesseur
            current = path[i]
            pred = path[i-1]
            
            # Vérifier si le nœud est dans L_M ET l'arc (predecessor, current) est dans A_M
            if current in L_M && (pred, current) in A_M
                push!(mtn_node, current)
            end
        end

        if round(Int, sum(y_val[j,ac] for j in mtn_node)) != 0.0
            # Trouver le nœud actuel de maintenance
            current_mtn = nothing
            for j in mtn_node
                if y_val[j,ac] >= 0.9
                    current_mtn = j
                    break
                end
            end
            
            current_idx = findfirst(==(current_mtn), path)
        
            # Prendre seulement les nœuds de maintenance AVANT current_mtn dans le chemin
            available_nodes = [j for j in mtn_node if findfirst(==(j), path) < current_idx]
            
            # Créer UNE solution pour CHAQUE nœud disponible
            for new_mtn in available_nodes
                #= if new_mtn == current_mtn
                    continue  # Sauter si c'est le même nœud
                end
                 =#
                y_new = copy(y_val)  # Copie de la solution actuelle
                y_new[current_mtn,ac] = 0  # Retirer l'ancienne maintenance
                y_new[new_mtn,ac] = 1      # Ajouter la nouvelle maintenance
                #= # Récupérer le chemin jusqu'à new_mtn
                new_mtn_idx = findfirst(==(new_mtn), path)
                path_to_new_mtn = path[1:new_mtn_idx] 
                 =#
                push!(y_list, y_new)  # Ajouter à la liste
                #println("Noeud de maintenance changé de $current_mtn à $new_mtn pour l'avion $ac")
            end
        end
    end
    return y_list
end

function spdual_anal(path, y_val, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
    result = (val_t1 = Dict(a => 0.0 for a in A_M), val_t2 = Dict(a => 0.0 for a in A_theta2) , val_s = Dict(j => 0.0 for j in L_M), 
                val_t3 = Dict(a => 0.0 for a in A_theta3), val_t4 = Dict(a => 0.0 for a in A_M), val_p = Dict(j => 0.0 for j in eq_nodes), 
                val_a1 = Dict(j => 0.0 for j in V), val_a2 = Dict(j => 0.0 for j in V))
    for i in 1:length(path)-1
        if (path[i], path[i+1]) in A_M && y_val[path[i+1]] == 1
            result.val_t1[(path[i], path[i+1])] = 1.0
        end
    end
end

function detect_i_succ_max(succ_A, V)
    max_succ = 0 
    i_max = nothing
    for i in setdiff(V, ["t"])
        temp = length(succ_A[i])
        println
        if temp >= max_succ
            max_succ = temp
            i_max = i
        end 
    end
    return i_max
end 

#= function x_neighborhood(x_val, i_max, pred_A, succ_A)
    temp = [i for i in pred_A[i_max] if x_val[(i, i_max)] >= 0.9]
    imax_pred_xval = isempty(temp) ? nothing : first(temp)
    println(succ_A[imax_pred_xval])
    for k in setdiff(succ_A[imax_pred_xval], ["s", "t"])
        if k != i_max
            println(k)
        end
    end 
end 
 =#

function x_neighborhood(x_val, fl_nodes, pred_A, succ_A, V)
    neighbors = []
    n_sample = min(5, length(fl_nodes))
    list_fl = rand(fl_nodes, n_sample)
    #list_fl = [detect_i_succ_max(succ_A, V)]
    for i in list_fl
        #println("i = $i")
        # Trouver le prédécesseur de i qui a un arc actif
        i_pred = [j for j in pred_A[i] if x_val[(j, i)] >= 0.9][1]
        isnothing(i_pred) && continue
        #println("i_pred = $i_pred")

        # Pour chaque successeur du prédécesseur (sauf s, t et i lui-même)
        for k in setdiff(succ_A[i_pred], ["s", "t", i])
            #println("\tk = $k")

            # Trouver le prédécesseur de k
            k_pred = [j for j in pred_A[k] if x_val[(j, k)] >= 0.9][1]
            isnothing(k_pred) && continue
            #println("\tk_pred = $k_pred")
            
            # Créer la solution voisine en échangeant i et k
            #println(pred_A[i])
            if k_pred in pred_A[i]
                x_new = copy(x_val)
                x_new[(i_pred, i)] = 0.0
                x_new[(i_pred, k)] = 1.0
                x_new[(k_pred, k)] = 0.0
                x_new[(k_pred, i)] = 1.0
                push!(neighbors, x_new)
                break
            end    
        end
    end
    
    return neighbors
end

function print_information(spd_result, instance, A_theta2, A_theta3, eq_nodes)
    A_M, L_M = instance["A_M"], instance["L_M"]
    d, max_flt = instance["d"], instance["maximum_flying_time"] 
    a_nodes, f, V = instance["a_nodes"], instance["initial_flying_time"], instance["V"]
    
    if spd_result.status == MOI.DUAL_INFEASIBLE
        println("\n VARAIBLE THETA1")
        for a in A_M
            i,j = a
            if spd_result.ray_t1[a] != 0    
                println("i = $i, j = $j : $(spd_result.ray_t1[a])")
            end
        end

        println("\n VARAIBLE THETA2")
        for a in A_theta2
            i,j = a
            if spd_result.ray_t2[a] != 0    
                println("i = $i, j = $j : $(spd_result.ray_t2[a])")
            end
        end 

        println("\n VARAIBLE SIGMA")
        for j in L_M
            if spd_result.ray_s[j] != 0    
                println("j = $j : $(spd_result.ray_s[j])")
            end 
        end

        println("\n VARAIBLE THETA3")
        for a in A_theta3
            i,j = a
            if spd_result.ray_t3[a] != 0    
                println("i = $i, j = $j : $(spd_result.ray_t3[a])")
            end
        end

        println("\n VARAIBLE THETA4")        
        for a in A_M
            i,j = a
            if spd_result.ray_t4[a] != 0    
                println("i = $i, j = $j : $(spd_result.ray_t4[a])")
            end
        end

        println("\n VARAIBLE PI")
        for k in eq_nodes
            if spd_result.ray_p[k] != 0    
                println("k = $k : $(spd_result.ray_p[k])")
            end
        end

        println("\n VARAIBLE ALPHA1")
        for j in Set(V)
            if spd_result.ray_a1[j] != 0    
                println("j = $j : $(spd_result.ray_a1[j])")
            end
        end

        println("\n VARAIBLE ALPHA2")
        for j in Set(V)
            if spd_result.ray_a2[j] != 0    
                println("j = $j : $(spd_result.ray_a2[j])")
            end
        end
    else
        println("\n VARAIBLE THETA1")
        for a in A_M
            i,j = a
            if spd_result.val_t1[a] != 0    
                println("i = $i, j = $j : $(spd_result.val_t1[a])")
            end
        end

        println("\n VARAIBLE THETA2")
        for a in A_theta2
            i,j = a
            if spd_result.val_t2[a] != 0    
                println("i = $i, j = $j : $(spd_result.val_t2[a])")
            end
        end 

        println("\n VARAIBLE SIGMA")
        for j in L_M
            if spd_result.val_s[j] != 0    
                println("j = $j : $(spd_result.val_s[j])")
            end 
        end

        println("\n VARAIBLE THETA3")
        for a in A_theta3
            i,j = a
            if spd_result.val_t3[a] != 0    
                println("i = $i, j = $j : $(spd_result.val_t3[a])")
            end
        end

        println("\n VARAIBLE THETA4")        
        for a in A_M
            i,j = a
            if spd_result.val_t4[a] != 0    
                println("i = $i, j = $j : $(spd_result.val_t4[a])")
            end
        end

        println("\n VARAIBLE PI")
        for k in eq_nodes
            if spd_result.val_p[k] != 0    
                println("k = $k : $(spd_result.val_p[k])")
            end
        end

        println("\n VARAIBLE ALPHA1")
        for j in Set(V)
            if spd_result.val_a1[j] != 0    
                println("j = $j : $(spd_result.val_a1[j])")
            end
        end

        println("\n VARAIBLE ALPHA2")
        for j in Set(V)
            if spd_result.val_a2[j] != 0    
                println("j = $j : $(spd_result.val_a2[j])")
            end
        end
    end 
end

function build_aircraft_path(start_aircraft::String, succ::Dict{String, String})
    """Construit le chemin complet d'un avion à partir du graphe des successeurs"""
    chemin = ["s", start_aircraft]
    current = start_aircraft
    
    while haskey(succ, current) && succ[current] != "t"
        current = succ[current]
        push!(chemin, current)
    end
    
    # Ajouter le nœud terminal si accessible
    haskey(succ, current) && push!(chemin, succ[current])
    
    return chemin
end 

function build_all_paths(x_val, A)
    """Construit les chemins de tous les avions à partir des variables x"""
    succ = Dict{String, String}()
    for (i, j) in Set(A)
        x_val[(i, j)] >= 0.9 && (succ[i] = j)
    end
    
    aircraft_paths = Dict{String, Vector{String}}()
    
    # Identifier les avions et construire leurs chemins
    aircraft_starts = [(i, j) for (i, j) in A if i == "s" && x_val[(i, j)] >= 0.9]
    
    for (_, avion) in aircraft_starts
        chemin = build_aircraft_path(avion, succ)
        aircraft_paths[avion] = chemin
    end
    
    return aircraft_paths
end

function find_aircraft_for_maintenance(maintenance_node::String, aircraft_paths::Dict{String, Vector{String}})
    """Trouve l'avion qui effectue la maintenance au nœud donné"""
    for (aircraft, path) in aircraft_paths
        maintenance_node in path && return aircraft
    end
    return "UNKNOWN"
end

function print_x_y(x_val, y_val, A, L_M, d, a_nodes)

    println("\nAFFICHAGE DE LA SOLUTION ")
    # Construire le graphe des successeurs en une seule passe
    succ = Dict{String, String}()
    for (i, j) in A
        x_val[(i, j)] >= 0.9 && (succ[i] = j)
    end
    
    # Extraire les chemins des avions et trouver le plus long
    aircraft_paths = Dict{String, Vector{String}}()
    longest_aircraft = ""
    longest_path = String[]
    longest_length = 0
    
    # Identifier les avions et construire leurs chemins
    aircraft_starts = [(i, j) for (i, j) in A if i == "s" && x_val[(i, j)] >= 0.9]
    
    println("✈️ Paths for each aircraft:")
    aircraft_mtn = Dict{String, Vector{String}}()
    for (_, avion) in aircraft_starts
        chemin = build_aircraft_path(avion, succ)
        aircraft_paths[avion] = chemin
        
        # Vérifier si c'est le chemin le plus long
        if length(chemin) > longest_length
            longest_length = calculate_path_distance(chemin, d)
            longest_path = chemin
            longest_aircraft = avion
        end
        cumul = 0
        path_with_cumul = String[]
        for i in 1:length(chemin)
            cumul += d[chemin[i]]
            push!(path_with_cumul, "$(chemin[i])_(d=$(d[chemin[i]]), Σ=$cumul)")
        end

        path_str = join(path_with_cumul, " ➡️  ")
        println("\n🛩️ Aircraft $avion path ($(length(chemin)-3) nodes): $path_str")
        println("   📊 Cumul total des d: $cumul")
        for j in L_M
            if y_val[j,avion] >= 0.9
                push!(get!(aircraft_mtn, avion, String[]), j)
            end
        end
    end
    println("")
    for ac in keys(aircraft_mtn)
        println("🛠️ Aircraft $ac performs maintenance at nodes: $(aircraft_mtn[ac])")
    end
    #println("🛠️ Maintenance at node $j by aircraft $avion")

    #= maintenance_flights = [j for j in L_M if sum(y_val[j,k] for k in a_nodes)  >= 0.9]
    maintenance_at_node = Dict{String, String}()
    
    for j in maintenance_flights
        aircraft = find_aircraft_for_maintenance(j, aircraft_paths)
        maintenance_at_node[j] = aircraft
        println("🛠️ Maintenance at node $j by aircraft $aircraft")
    end =#
end

#= function find_cover_set(A, succ_A, x_val, a_nodes, d, max_flt, L_M)
    succ = Dict{String, String}()
    for (i, j) in Set(A)
        x_val[(i, j)] >= 0.9 && (succ[i] = j)
    end
    all_cover_set = Vector{Tuple{Vector{String}, Int}}()
    for aircraft in a_nodes
        cover_set = ["s", aircraft]
        current = aircraft
        temp = d[aircraft]
        len_cvs = 2
        while true
            next = get(succ, current, nothing)
            next == nothing || next == "t" && break
            temp += d[next]
            push!(cover_set, next)
            current = next
            len_cvs += 1
            #println("\nNEXT : $next \tTEMP = $temp")
            temp > max_flt && break
        end
        #= if sum(y_val[j] for j in intersect(L_M, cover_set)) == 0
            push!(all_cover_set, (copy(cover_set), len_cvs))
        end =#
        if temp > max_flt        
            #= println("\nyesTEMP FINAL = $temp")
            println("LEN = $len_cvs")
            println("COVER SET = $cover_set")
            println(all_cover_set)
            =#
            key_node = cover_set[end-1]
            temp -= d[cover_set[end]]
            pop!(cover_set)
            
            for j in get(succ_A, key_node, [])
                if temp + d[j] > max_flt
                    push!(cover_set, j)
                    #= println("\nTEMP FINAL = $temp")
                    println("LEN = $len_cvs")
                    println("COVER SET = $cover_set")
                     =#
                    push!(all_cover_set, (copy(cover_set), len_cvs))
                    #= if sum(y_val[j] for j in intersect(L_M, cover_set)) == 0
                        push!(all_cover_set, (copy(cover_set), len_cvs))
                    end 
                     =#
                    pop!(cover_set)
                end
            end
        end 
        
    end    
    return all_cover_set
end
 =#

function find_cover_set(A, succ_A, x_val, y_val, a_nodes, d, max_flt, L_M)
    succ = Dict{String, String}()
    sizehint!(succ, length(A))
    @inbounds for (i, j) in A
        if x_val[(i, j)] >= 0.9
            succ[i] = j
        end
    end
    # 2. Pré-calculer L_M comme Set pour intersections rapides
    L_M_set = Set(L_M)
    
    # 3. Pré-allouer le vecteur de résultats avec estimation de taille
    all_cover_set = Vector{Tuple{Vector{String}, Int, String}}()
    sizehint!(all_cover_set, length(a_nodes) * 2)
    
    # 4. Réutiliser le buffer cover_set au lieu de créer/détruire
    cover_set = Vector{String}()
    sizehint!(cover_set, 20)  # Estimation taille typique
    
    @inbounds for aircraft in a_nodes
        # Réinitialiser au lieu de recréer
        empty!(cover_set)
        push!(cover_set, "s", aircraft)
        
        current = aircraft
        temp = d[aircraft]
        len_cvs = 2
        
        # 5. Suivre le chemin avec get optimisé
        while true
            next = get(succ, current, nothing)
            
            # Break conditions combinées
            if isnothing(next) || next == "t"
                break
            end
            
            temp += d[next]
            push!(cover_set, next)
            len_cvs += 1
            current = next
            
            # Sortir dès qu'on dépasse
            if temp > max_flt
                if sum(y_val[j,aircraft] for j in intersect(L_M_set, cover_set); init=0.0) >= 0.9
                    break
                #= else
                    push!(all_cover_set, (copy(cover_set), len_cvs))  =#
                end
                 # 6. Traitement de la violation
                key_node = cover_set[end-1]
                temp -= d[cover_set[end]]
                pop!(cover_set)
                
                # 7. Obtenir successeurs une seule fois
                successors = get(succ_A, key_node, nothing)
                for j in successors
                    if temp + d[j] > max_flt
                        push!(cover_set, j)
                        # Copier seulement si nécessaire
                        push!(all_cover_set, (copy(cover_set), len_cvs, aircraft))
                        
                        pop!(cover_set)
                    end
                end 
                break
            end
        end
    end
    
    return all_cover_set
end

#= function find_cover_set_2(fl_nodes, a_nodes, d, max_flt, pred_A, succ_A)
    nodes_tries = sort(union(fl_nodes, a_nodes), by = x -> d[x], rev = true)
    #= cv_set = []
    somme = 0
    len_cvs = 0
     =#
    all_cv_set = []
    for node in a_nodes
        cv_set = [node]
        somme = d[node]
        len_cvs = 1
        while somme < max_flt
            if cv_set[end] == "t"
                break
            end
            next_nodes = sort(succ_A[cv_set[end]], by = x -> d[x], rev = true)
            temp_node = next_nodes[1]
            if !(temp_node in cv_set)
                push!(cv_set, temp_node)
                len_cvs += 1
                somme += d[temp_node]
            end
        end 
        push!(all_cv_set, (copy(cv_set), len_cvs, somme))
    end
    return all_cv_set
end 
=#
#= function find_cover_set_2(fl_nodes, a_nodes, d, max_flt, pred_A, succ_A, max_per_node=10)
    all_cv_set = []
    
    # Fonction récursive pour explorer les chemins
    function explore_paths(cv_set, somme, len_cvs, visited, results)
        # Si on a déjà assez de résultats pour ce nœud, arrêter
        if length(results) >= max_per_node
            return
        end
        
        # Si on dépasse max_flt, on enregistre cet ensemble
        if somme > max_flt
            push!(results, (copy(cv_set), len_cvs, somme))
            return
        end
        
        # Si on atteint "t", on arrête ce chemin
        if !isempty(cv_set) && cv_set[end] == "t"
            return
        end
        
        # Obtenir les successeurs du dernier nœud
        current_node = cv_set[end]
        
        if haskey(succ_A, current_node)
            next_nodes = succ_A[current_node]
            
            # Trier les successeurs par valeur décroissante (heuristique gloutonne)
            sorted_next = sort(collect(next_nodes), by = x -> d[x], rev = true)
            
            # Explorer les successeurs
            for next_node in sorted_next
                # Arrêter si on a assez de résultats
                if length(results) >= max_per_node
                    break
                end
                
                if !(next_node in visited)  # Éviter les cycles
                    push!(cv_set, next_node)
                    push!(visited, next_node)
                    
                    explore_paths(cv_set, somme + d[next_node], len_cvs + 1, visited, results)
                    
                    pop!(cv_set)
                    pop!(visited)
                end
            end
        end
    end
    
    # Démarrer l'exploration depuis chaque nœud de départ dans a_nodes
    for node in a_nodes
        results_per_node = []
        cv_set = [node]
        visited = Set([node])
        somme = d[node]
        len_cvs = 1
        
        explore_paths(cv_set, somme, len_cvs, visited, results_per_node)
        
        # Ajouter les résultats de ce nœud au résultat global
        append!(all_cv_set, results_per_node)
    end
    
    return all_cv_set
end
 =#
function find_cover_set_min(a_nodes, d, max_flt, succ_A)
    all_cv_set = []
    
    for start_node in a_nodes
        min_length = Inf
        node_results = []
        memo = Dict()  # Mémorisation
        
        function explore_paths(cv_set, somme, len_cvs, visited)
            # Clé pour mémorisation
            state_key = (cv_set[end], somme, len_cvs)
            if haskey(memo, state_key)
                return
            end
            memo[state_key] = true
            
            if somme > max_flt
                if len_cvs < min_length
                    min_length = len_cvs
                    empty!(node_results)
                    push!(node_results, (copy(cv_set), len_cvs, somme))
                elseif len_cvs == min_length
                    push!(node_results, (copy(cv_set), len_cvs, somme))
                end
                return
            end
            
            if len_cvs >= min_length
                return
            end
            
            current_node = cv_set[end]
            if current_node == "t"
                return
            end
            
            if haskey(succ_A, current_node)
                for next_node in succ_A[current_node]
                    if !(next_node in visited)
                        push!(cv_set, next_node)
                        push!(visited, next_node)
                        
                        explore_paths(cv_set, somme + d[next_node], len_cvs + 1, visited)
                        
                        pop!(cv_set)
                        pop!(visited)
                    end
                end
            end
        end
        
        cv_set = [start_node]
        visited = Set([start_node])
        explore_paths(cv_set, d[start_node], 1, visited)
        
        append!(all_cv_set, node_results)
    end
    
    return all_cv_set
end

function find_cover_set_2(a_nodes, d, max_flt, succ_A, nmax = 50)
    all_cv_set = []
    
    for start_node in a_nodes
        node_results = []
        count = 0
        memo = Dict()  # Mémorisation
        
        function explore_paths(cv_set, somme, len_cvs, visited)
            # Si on a déjà 100 chemins, arrêter
            if count >= 100
                return true
            end
            
            # Clé pour mémorisation (nœud actuel + somme)
            state_key = (cv_set[end], somme)
            if haskey(memo, state_key)
                return false
            end
            memo[state_key] = true
            
            # Si on dépasse max_flt, on enregistre
            if somme > max_flt
                push!(node_results, (copy(cv_set), len_cvs, somme))
                count += 1
                
                if count >= nmax
                    return true
                end
                return false  # Arrêter cette branche
            end
            
            # Obtenir les successeurs du dernier nœud
            current_node = cv_set[end]
            next_nodes = get(succ_A, current_node, nothing)
            
            if next_nodes !== nothing
                # Explorer tous les successeurs possibles
                for next_node in next_nodes
                    if next_node ∉ visited  # Éviter les cycles
                        push!(cv_set, next_node)
                        push!(visited, next_node)
                        
                        # Explorer et vérifier si on doit arrêter
                        should_stop = explore_paths(cv_set, somme + d[next_node], len_cvs + 1, visited)
                        
                        pop!(cv_set)
                        delete!(visited, next_node)
                        
                        # Si on a atteint 100 chemins, arrêter l'exploration
                        if should_stop
                            return true
                        end
                    end
                end
            end
            return false
        end
        
        cv_set = [start_node]
        visited = Set([start_node])
        explore_paths(cv_set, d[start_node], 1, visited)
        
        append!(all_cv_set, node_results)
    end
    
    return all_cv_set
end
#= function find_cover_set_2(a_nodes, d, max_flt, succ_A)
    all_cv_set = []
    #Cherche tous les ensembles de couverture possibles en dépassant max_flt
    # Fonction récursive pour explorer les chemins
    function explore_paths(cv_set, somme, len_cvs, visited)
        # Si on dépasse max_flt, on enregistre et on arrête cette branche
        if somme > max_flt
            push!(all_cv_set, (copy(cv_set), len_cvs, somme))
            return  # Arrêter ici, pas besoin de continuer
        end
        
        # Obtenir les successeurs du dernier nœud
        current_node = cv_set[end]
        
        if haskey(succ_A, current_node)
            next_nodes = succ_A[current_node]
            
            # Explorer tous les successeurs possibles
            for next_node in next_nodes
                if !(next_node in visited)  # Éviter les cycles
                    push!(cv_set, next_node)
                    push!(visited, next_node)
                    
                    explore_paths(cv_set, somme + d[next_node], len_cvs + 1, visited)
                    
                    pop!(cv_set)
                    pop!(visited)
                end
            end
        end
    end
    
    # Démarrer l'exploration depuis chaque nœud de départ dans a_nodes
    for node in a_nodes
        cv_set = [node]
        visited = Set([node])
        somme = d[node]
        len_cvs = 1
        
        explore_paths(cv_set, somme, len_cvs, visited)
    end
    
    return all_cv_set
end
 =#
#= function find_cover_set_2(fl_nodes, d, max_flt)
    nodes_tries = sort(fl_nodes, by = x -> d[x], rev = true)
    cv_set = []
    somme = 0
    all_cv_set = []
    len_cvs = 0
    
    for node in nodes_tries
        push!(cv_set, node)
        len_cvs += 1
        somme += d[node]
        
        if somme > max_flt
            push!(all_cv_set, (copy(cv_set), len_cvs, somme))  # ← Ajout de copy()
            pop!(cv_set)
            len_cvs -= 1
            somme -= d[node]
        end
    end
    
    return all_cv_set
end
 =#
function find_mtn_set(all_cover_set, A_M, d, max_flt)
    #all_mtn_set = Vector{NamedTuple{(:set_node, :len_set, :bnd_val), Tuple{Vector{String}, Int, Int}}}()
    best_mtn_set = (set_node = [], len_set = 0, bnd_val = 0)
    best_bnd_val = max_flt
    for (cover_set, len_cvs) in all_cover_set
        mtn_set = copy(cover_set)
        len_mtns = copy(len_cvs)
        end_node = mtn_set[end]
        while !((mtn_set[end-1], end_node) in A_M) && len_mtns > 2
            pop!(mtn_set)
            len_mtns -= 1
            end_node = mtn_set[end]
        end
        bnd_val = max_flt -  compute_fl_time(mtn_set, d) + d[end_node]
        #println("\nMAINTENANCE SET = $mtn_set")
        #println("BOUND = $bnd_val")
        #push!(all_mtn_set, (set_node = copy(mtn_set), len_set = len_mtns, bnd_val = bnd_val))
        if bnd_val < best_bnd_val
            best_bnd_val = bnd_val
            best_mtn_set = (set_node = copy(mtn_set), len_set = len_mtns, bnd_val = bnd_val)
        end
    end
    return best_mtn_set
end 

function compute_fl_time(set_nodes, d)
    total_time = 0
    for node in set_nodes
        total_time += d[node]
    end
    return total_time
end

