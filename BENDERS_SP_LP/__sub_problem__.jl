include("functions_for_benders.jl")

function build_spdual(env, instance_data)
    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data
    other_data = instance_data.other_data

    L_M = node_sets.mtn_nodes
    A_M = arc_sets.arcs_M
    A_theta2 = other_data["A_theta2"]
    A_theta3 = other_data["A_theta3"]
    eq_nodes = other_data["eq_nodes"]
    f = fl_data.init_flying_time
    V = graph.nodes
    d = fl_data.d

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

function solve_spdual(x_val, y_val, instance_data, spdual_model)
    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data
    other_data = instance_data.other_data
    
    L_M = node_sets.mtn_nodes
    A_M = arc_sets.arcs_M
    f = fl_data.init_flying_time
    d = fl_data.d
    a_nodes = node_sets.ac_nodes
    V = graph.nodes
    max_flt = instance_data.max_flying_time

    A_theta2 = other_data["A_theta2"]
    A_theta3 = other_data["A_theta3"]
    
    # Récupérer les variables du modèle dual
    theta1, theta2, sigma = spdual_model[:theta1], spdual_model[:theta2], spdual_model[:sigma]
    theta3, theta4, pi1 = spdual_model[:theta3], spdual_model[:theta4], spdual_model[:pi1]
    alpha1, alpha2 = spdual_model[:alpha1], spdual_model[:alpha2]
    
    # Mettre à jour l'objectif avec les valeurs courantes de x et y
    @objective(spdual_model, Max,   
        sum((max_flt*x_val[(i,j)] - (max_flt - d[i])*(1 - y_val[j])) * theta1[(i,j)] for (i,j) in A_M) + 
        sum((d[j] + (max_flt - d[i] - d[j])*(1 - x_val[(i,j)])) * theta2[(i,j)] for (i,j) in A_theta2) +
        sum((max_flt - (max_flt - d[j])*y_val[j]) * sigma[j] for j in L_M) +
        sum((d[j] - max_flt*(1 - x_val[(i,j)])) * theta3[(i,j)] for (i,j) in A_theta3) +
        sum((d[j] - max_flt*(1 - x_val[(i,j)]) - max_flt*y_val[j]) * theta4[(i,j)] for (i,j) in A_M) + 
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

function fix_y_decision(aircrafts_paths, instance_data)
    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data

    A = graph.arcs
    d = fl_data.d
    max_flt = instance_data.max_flying_time
    H_aircraft = node_sets.H_aircraft
    L_M = node_sets.mtn_nodes
    A_M = arc_sets.arcs_M

    y_val = Dict(j => 0.0 for j in L_M)
    mtn_nodes = Dict(aircraft => [] for aircraft in H_aircraft)
    status = "FAISABLE"
    for aircraft in H_aircraft
        path = aircraft_paths[aircraft]
        len_p = length(path)
        summ = d[path[1]]
        len_mtn_p = 0
        for i in 2:len_p
            if (path[i-1], path[i]) in A_M
                push!(mtn_nodes[aircraft], path[i]) 
            elseif summ + d[path[i]] > max_flt
                break
            end 
            summ += d[path[i]]
        end 
    end
    for aircraft in H_aircraft
        y
    end 

end

function find_subpath_mtn(x_val, result, instance_data)
    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data

    A = graph.arcs
    d = fl_data.d
    max_flt = instance_data.max_flying_time
    H_aircraft = node_sets.H_aircraft
    L_M = node_sets.mtn_nodes
    A_M = arc_sets.arcs_M


    y_val = result.y
    aircraft_paths = build_all_paths(x_val, A)
    mtn_sub_path = Dict{String, Tuple{Vector{String}, Int}}()
    infeas_sub_path = Dict{String, Tuple{Vector{String}, Int}}()
    for aircraft in H_aircraft
        path = aircraft_paths[aircraft]
        len_p = length(path)
    
        if result.status == MOI.OPTIMAL
            for i in 1:len_p
                if path[i] in L_M && y_val[path[i]] >= 0.9
                    mtn_sub_path[aircraft] = (path[1:i], i)
                    break
                end 
            end
        elseif result.status == MOI.INFEASIBLE
            mtn_poss = false
            summ = d[path[1]]
            for i in 2:len_p
                if (path[i-1], path[i]) in A_M
                    #= println("\n----------UN ARC DE MAINTENANCE EST TROUVÉ ----------\n")
                    println(path[i-1], "  ",  path[i])
                     =#
                    break  
                elseif summ + d[path[i]] > max_flt
                    infeas_sub_path[aircraft] = (path[1:i], i)
                    break
                end 
                summ += d[path[i]]
            end
        end 
    end 

    if result.status == MOI.OPTIMAL
        return mtn_sub_path
    elseif result.status == MOI.INFEASIBLE
        return infeas_sub_path
    end 

end


