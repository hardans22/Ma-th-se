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
        sum((max_flt*sum(x_val[(i,j,k)] for k in a_nodes) - (max_flt - d[i])*(1 - sum(y_val[j,k] for k in a_nodes))) * theta1[(i,j)] for (i,j) in A_M) + 
        sum((d[j] + (max_flt - d[i] - d[j])*(1 - sum(x_val[(i,j,k)] for k in a_nodes))) * theta2[(i,j)] for (i,j) in A_theta2) +
        sum((max_flt - (max_flt - d[j])*sum(y_val[j,k] for k in a_nodes)) * sigma[j] for j in L_M) +
        sum((d[j] - max_flt*(1 - sum(x_val[(i,j,k)] for k in a_nodes))) * theta3[(i,j)] for (i,j) in A_theta3) + 
        sum((d[j] - max_flt * (1 - sum(x_val[(i,j,k)] for k in a_nodes)) - max_flt*sum(y_val[j,k] for k in a_nodes))*theta4[(i,j)] for (i,j) in A_M) + 
        sum(f[k] * pi1[k] for k in a_nodes) + sum(d[j] * alpha1[j] for j in V) + sum(max_flt * alpha2[j] for j in V))

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

function solve_aux_spdual(x_val, y_val, x_core, y_core, instance, aux_spdual, sp_obj)
    A_M, L_M, A_theta2, A_theta3 = instance["A_M"], instance["L_M"], instance["A_theta2"], instance["A_theta3"]
    d, max_flt, eq_nodes = instance["d"], instance["maximum_flying_time"], instance["eq_nodes"] 
    a_nodes, f, V = instance["a_nodes"], instance["initial_flying_time"], instance["V"]
    
    # Récupérer les variables du modèle dual
    theta1, theta2, sigma = aux_spdual[:theta1], aux_spdual[:theta2], aux_spdual[:sigma]
    theta3, theta4, pi1 = aux_spdual[:theta3], aux_spdual[:theta4], aux_spdual[:pi1]
    alpha1, alpha2 = aux_spdual[:alpha1], aux_spdual[:alpha2]
    
    # Mettre à jour l'objectif avec les valeurs courantes de x et y
    @objective(aux_spdual, Max,   
        sum((max_flt*sum(x_core[(i,j,k)] for k in a_nodes) - (max_flt - d[i])*(1 - sum(y_core[j,k] for k in a_nodes))) * theta1[(i,j)] for (i,j) in A_M) + 
        sum((d[j] + (max_flt - d[i] - d[j])*(1 - sum(x_core[(i,j,k)] for k in a_nodes))) * theta2[(i,j)] for (i,j) in A_theta2) +
        sum((max_flt - (max_flt - d[j])*sum(y_core[j,k] for k in a_nodes)) * sigma[j] for j in L_M) +
        sum((d[j] - max_flt*(1 - sum(x_core[(i,j,k)] for k in a_nodes))) * theta3[(i,j)] for (i,j) in A_theta3) +
        sum((d[j] - max_flt*(1 - sum(x_core[(i,j,k)] for k in a_nodes)) - max_flt*sum(y_core[j,k] for k in a_nodes)) * theta4[(i,j)] for (i,j) in A_M) + 
        sum(f[k] * pi1[k] for k in a_nodes) + sum(d[j] * alpha1[j] for j in V) +  sum(max_flt * alpha2[j] for j in V))

    if aux_spdual.ext[:c5_ref] !== nothing
        delete(aux_spdual, aux_spdual.ext[:c5_ref])
    end

    aux_spdual.ext[:c5_ref] = @constraint(aux_spdual,
        sum((max_flt*sum(x_val[(i,j,k)] for k in a_nodes) - (max_flt - d[i])*(1 - sum(y_val[j,k] for k in a_nodes))) * theta1[(i,j)] for (i,j) in A_M) + 
        sum((d[j] + (max_flt - d[i] - d[j])*(1 - sum(x_val[(i,j,k)] for k in a_nodes))) * theta2[(i,j)] for (i,j) in A_theta2) +
        sum((max_flt - (max_flt - d[j])*sum(y_val[j,k] for k in a_nodes)) * sigma[j] for j in L_M) +
        sum((d[j] - max_flt*(1 - sum(x_val[(i,j,k)] for k in a_nodes))) * theta3[(i,j)] for (i,j) in A_theta3) +
        sum((d[j] - max_flt*(1 - sum(x_val[(i,j,k)] for k in a_nodes)) - max_flt*sum(y_val[j,k] for k in a_nodes)) * theta4[(i,j)] for (i,j) in A_M) + 
        sum(f[k] * pi1[k] for k in a_nodes) + sum(d[j] * alpha1[j] for j in V) + 
        sum(max_flt * alpha2[j] for j in V) == sp_obj)

    optimize!(aux_spdual)
    spdual_status = termination_status(aux_spdual)
    #println("Aux SP Dual status: ", spdual_status)
    
    if spdual_status == MOI.OPTIMAL
        #sprintln("sp_obj: ", sp_obj)
        return (status = spdual_status, obj = objective_value(aux_spdual), 
                val_t1 = value.(theta1), val_t2 = value.(theta2), val_s = value.(sigma), 
                val_t3 = value.(theta3), val_t4 = value.(theta4), val_p = value.(pi1), 
                val_a1 = value.(alpha1), val_a2 = value.(alpha2), 
                time = round(solve_time(aux_spdual), digits = 6))
    else
        println("\nAux SP Dual status: ", spdual_status)
        println("sp_obj: ", sp_obj)
        return (status = spdual_status, obj = nothing)
    end
end 

function master_LP(env, A, A_T, L_M, V_wt_st, a_nodes, succ_A, pred_A, pred_A_M, acritik_set, ac_critique)
    core_model = Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(core_model, "OutputFlag", 0)
    set_optimizer_attribute(core_model, "Threads", nbr_thread)
    # Variables
    @variable(core_model, x[a in A, k in a_nodes], Bin)
    @variable(core_model, y[j in L_M, k in a_nodes], Bin)
    @variable(core_model, r >= 0)

    NH_aircraft = setdiff(a_nodes, acritik_set)
    @constraint(core_model, mp_c1[k in a_nodes], x[("s",k),k] == 1)
    @constraint(core_model, mp_c2[k in a_nodes], sum(x[a,k] for a in A_T) == 1)
    @constraint(core_model, mp_c3[i in V_wt_st], sum(x[(i,j),k] for j in get(succ_A, i, []), k in a_nodes) == 1)
    @constraint(core_model, mp_c4[i in V_wt_st, k in a_nodes], sum(x[(j,i),k] for j in get(pred_A, i, [])) == sum(x[(i,j),k] for j in get(succ_A, i, [])))
    @constraint(core_model, mp_c6[j in L_M, k in NH_aircraft], sum(y[j,k] for j in L_M) == 0)

    #= @constraint(core_model, mp_c5[j in L_M, k in a_nodes], y[j,k] + r <= sum(x[(i,j),k] for i in get(pred_A_M, j, [])))
    @constraint(core_model, sum(y[j,k] for j in L_M, k in a_nodes) + r <= ac_critique)
    @constraint(core_model, mp_c7[k in acritik_set], sum(y[j,k] for j in L_M) + r <= 1)
     =##@constraint(core_model, mp_c8[j in L_M], sum(y[j,k] for k in a_nodes) <= 1)

     
    @objective(core_model, Min, r)

    optimize!(core_model)

    x_val = Dict((i,j,k) => value(x[(i,j),k]) for (i,j) in A, k in a_nodes)
    y_val = Dict((j,k) => value(y[j,k]) for j in L_M, k in a_nodes)
    r_val = value(r)
    for k in acritik_set
        for j in L_M
            if sum(x_val[(i,j,k)] for i in get(pred_A_M, j, [])) != 0.0
                y_val[(j,k)] = sum(x_val[(i,j,k)] for i in get(pred_A_M, j, []))/(length(L_M))
            end
        end
    end 
    return x_val, y_val
end

function find_interior_point1(x_val, y_val, A, L_M, a_nodes)
    alpha = 0.1  # 20% vers le centre, 80% du point extrême
    
    #= x_core = Dict((i,j,k) => (1-alpha) * x_val[(i,j,k)] + alpha * 0.5 
                 for (i,j) in A, k in a_nodes) =#
    
    y_core = Dict((j,k) => (1-alpha) * y_val[(j,k)] + alpha * 0.5 
                 for j in L_M, k in a_nodes)
    
    return x_val, y_core
end

function find_interior_point2(x_val, y_val, A, L_M, a_nodes) 
                              
    alpha=0.2
    """
    Combine le point extrême avec un point "central" artificiel
    x_interior = (1-α)·x_extreme + α·x_center
    """
    
    x_core = Dict()
    for (i,j) in A, k in a_nodes
        # Point central = 0.5 pour toutes les variables continues
        x_center = 0.5
        x_core[(i,j,k)] = (1 - alpha) * x_val[(i,j,k)] + alpha * x_center
    end
    
    y_core = Dict()
    for j in L_M, k in a_nodes
        y_center = 0.5
        y_core[(j,k)] = (1 - alpha) * y_val[(j,k)] + alpha * y_center
    end
    
    return x_core, y_core
end

function MP_is_feasible(x_val, y_val, master_model, a_nodes, A, L_M)
    mp_copy, ref_map = copy_model(master_model)
    set_optimizer(mp_copy, Gurobi.Optimizer)
    set_optimizer_attribute(mp_copy, "OutputFlag", 0)
    var_names = [:x, :y]
    for var_name in var_names
        mp_copy[var_name] = ref_map[master_model[var_name]]
    end
    for k in a_nodes 
        for (i,j) in A
            fix(mp_copy[:x][(i,j),k], x_val[(i,j,k)])
        end 
        for j in L_M
            fix(mp_copy[:y][j,k], y_val[(j,k)])
        end
    end 
    optimize!(mp_copy)
    mp_status = termination_status(mp_copy)
    println("MP copy status: ", mp_status)
    if mp_status != MOI.OPTIMAL
        return false
    end
    return true
end

function sp_analytic(ac_paths, ac_mtn, d, L_M, max_flt)
    run_time = @elapsed begin
        solution = Dict{String, Vector{String}}()
        s_rho = Dict(j => 0 for j in L_M)
        #println(s_rho)
        su = Dict()
        obj_val = 0.0 
        for (ac, path) in ac_paths
            su[ac] = d[ac]
            mtn_node = ac_mtn[ac]
            for i in 2:length(path)-1
                if path[i] == "s" || path[i] == "t" || path[i+1] == "t"
                    continue
                end
                if path[i+1] in L_M 
                    if path[i+1] == mtn_node
                        s_rho[path[i+1]] = max_flt - su[path[i]]
                        su[path[i+1]] = d[path[i+1]]
                        solution[ac] = path[1:i+1]
                    else
                        s_rho[path[i+1]] = 0
                        su[path[i+1]] = su[path[i]] + d[path[i+1]] 
                    end
                    if s_rho[path[i+1]] < 0
                        solution[ac] = path[1:i+1]
                        #println("\n\n IIIIIIIIIIINFEASIBLEEEEEEEEEEEEEEEE")
                    end
                    obj_val += s_rho[path[i+1]]
                else 
                    su[path[i+1]] = su[path[i]] + d[path[i+1]] 
                end
                
            end
            #su["t"] = su[path[end-1]]
        end 
    end 
    return (su = su, s_rho = s_rho, obj_val = obj_val, time = run_time, sol = solution)
end 