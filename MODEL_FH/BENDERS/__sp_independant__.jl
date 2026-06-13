mutable struct SPSets
    rdc_V::Vector{String}
    rdc_L_M::Vector{String}
    rdc_A::Vector{Tuple{String,String}}
    rdc_A_M::Vector{Tuple{String,String}}
    rdc_A_M_bar::Vector{Tuple{String,String}}
end 

function set_for_SP(aircraft_paths, instance_data)
    graph       = instance_data.graph
    node_sets   = graph.node_sets
    arc_sets    = graph.arc_sets
    fl_data     = instance_data.fl_data

    L_M = node_sets.mtn_nodes
    V   = graph.nodes
    A_M = arc_sets.arcs_M
    a_nodes     = node_sets.ac_nodes

    sets_for_SP = Dict()

    for ac in a_nodes
        path = aircraft_paths[ac].path
        len_p = aircraft_paths[ac].len_p
        rdc_V = copy(path)
        rdc_A = [(path[i], path[i+1]) for i in 1:len_p-1]
        rdc_L_M = filter(i -> i in L_M, rdc_V)
        rdc_A_M = filter(i -> i in A_M, rdc_A)
        rdc_A_M_bar = setdiff(rdc_A, rdc_A_M)
        sets_for_SP[ac] = SPSets(rdc_V, rdc_L_M, rdc_A, rdc_A_M, rdc_A_M_bar)
    end
    return sets_for_SP 
end 

function build_sp_ind(env, aircraft, set_for_SP, instance_data)
    fl_data     = instance_data.fl_data

    max_flt = instance_data.max_flying_time
    f       = fl_data.init_flying_time
    d       = fl_data.d

    rdc_V       = set_for_SP.rdc_V
    rdc_A       = set_for_SP.rdc_A
    rdc_A_M     = set_for_SP.rdc_A_M
    rdc_A_M_bar = set_for_SP.rdc_A_M_bar

    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(rdc_A_M, pred_A_M, succ_A_M)
    
    sp_model = Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(sp_model, "OutputFlag", 0)
    set_optimizer_attribute(sp_model, "InfUnbdInfo", 1)
    set_optimizer_attribute(sp_model, "DualReductions", 0)

    # ===================Decision variables ===================
    sp_model[:x_copy] = @variable(sp_model, x_copy[k in rdc_A])                         #Copy of x
    sp_model[:y] = @variable(sp_model, y[j in rdc_L_M], Bin)        
    sp_model[:u] = @variable(sp_model, d[j] <= u[j in rdc_V] <= max_flt)                #Accumulatad flying time at node j  
    sp_model[:rho] = @variable(sp_model, rho[j in rdc_L_M] >= 0)                        #Remaining flying time at node j

    # =========== Objective ===========
    @objective(sp_model, Min,  sum(rho[j] for j in rdc_L_M))

    # ========  Flying time constraints ========

    @constraint(sp_model, sum(y[j] for j in rdc_L_M) <= 1)
    @constraint(sp_model, c3[j in rdc_L_M], y[j] <= sum(x_copy[(i,j)] for i in get(pred_rdc_A_M, j, [])))
    @constraint(sp_model, c4[(i,j) in rdc_A_M], rho[j] >= max_flt*x_copy[(i, j)] - u[i] - (max_flt - d[i])*(1 - y[j]))
    @constraint(sp_model, c5[(i,j) in rdc_A; j != "t"], u[j] <= u[i] + d[j] + (max_flt - d[i] - d[j])*(1 - x_copy[(i,j)]))
    @constraint(sp_model, c6[j in rdc_L_M], u[j] <= max_flt -(max_flt - d[j])*y[j])
    @constraint(sp_model, c7[(i,j) in rdc_A_M_bar; j != "t"], u[j] >= u[i] + d[j] - max_flt*(1 - x_copy[(i,j)]))
    @constraint(sp_model, c8[(i,j) in rdc_A_M], u[j] >= u[i] + d[j] - max_flt*(1 - x_copy[(i,j)]) - max_flt*y[j])
    @constraint(sp_model, u["s"] == 0)
    @constraint(sp_model, u[k] == f[k])

    return sp_model
end

function solve_sp_ind(sp_model, x_val, rdc_L_M, rdc_V)
    # Fixer x_copy
    x_copy = sp_model[:x_copy]
    for arc in keys(x_val)
        fix(x_copy[arc], abs(x_val[arc]); force = true)
    end

    optimize!(sp_model)
    sp_status = termination_status(sp_model)
    if sp_status == MOI.OPTIMAL
        u_val = Dict(j => value(sp_model[:u][j]) for j in rdc_V)
        y_val = Dict(j => value(sp_model[:y][j]) for j in rdc_L_M)
        rho_val = Dict(j => value(sp_model[:rho][j]) for j in rdc_L_M)
        return (status = "OPTIMAL", obj = round(objective_value(sp_model)), 
                time = round(solve_time(sp_model), digits = 6), y = y_val, 
                u = u_val, rho = rho_val)
    else sp_status == MOI.INFEASIBLE
        return (status = "INFEASIBLE", obj = Inf, 
                y = Dict(j => 0.0 for j in rdc_L_M), 
                time = round(solve_time(sp_model), digits = 6))
    end
end 

function solve_sp_H_ind(aircraft, set_for_SP, ac_path, instance_data)
    fl_data     = instance_data.fl_data

    max_flt = instance_data.max_flying_time
    f       = fl_data.init_flying_time
    d       = fl_data.d

    rdc_V       = set_for_SP.rdc_V
    rdc_L_M     = set_for_SP.rdc_L_M
    rdc_A_M     = set_for_SP.rdc_A_M

    A_M_set      = Set(rdc_A_M)                  
    infeasible = false

    time = @elapsed begin
        u   = Dict{String, Float64}(j => 0.0 for j in rdc_V)
        y   = Dict{String, Int}(j => 0 for j in rdc_L_M)
        rho = Dict{String, Float64}(j => 0.0 for j in rdc_L_M)
        obj = 0.0
        status = "OPTIMAL"
        path = ac_path.[aircraft].path
        len_path = ac_path.len_path
        u[aircraft] = f[aircraft]
        mtn_node = nothing
        pred_mtn_node = nothing
        mtn_node_idx = nothing
        idx = 2 
        while idx < len_path-1
            node_cur  = path[idx]
            node_next = path[idx+1]

            if (node_cur, node_next) in A_M_set
                mtn_node      = node_next
                pred_mtn_node = node_cur
                mtn_node_idx  = idx + 1
            end
            if u[node_cur] + d[node_next] > max_flt
                if mtn_node !== nothing
                    y[mtn_node]   = 1
                    rho[mtn_node] = max_flt - u[pred_mtn_node]
                    u[mtn_node]   = d[mtn_node]
                    idx           = mtn_node_idx
                    obj          += rho[mtn_node]
                    continue
                else
                    status = "INFEASIBLE"
                    infeasible = true
                    break
                end
            end

            u[node_next] = u[node_cur] + d[node_next]
            idx += 1
        end  
        u["t"] = 0.0
    end 
    if infeasible
        return (status = status, obj = Inf, time = round(time, digits = 6))
    end 
    return (status = status, obj = obj, y = y, u = u, 
            rho = rho, time = round(time, digits = 6))
end 


function irreductible_path_ind(ac_path, set_for_SP, result, instance_data)
    graph      = instance_data.graph
    fl_data    = instance_data.fl_data
    d          = fl_data.d
    max_flt    = instance_data.max_flying_time

    L_M_set = Set(set_for_SP.rdc_L_M)
    A_M_set = Set(set_for_SP.rdc_A_M)

    if result.status == "OPTIMAL"
        y_val = result.y
        mtn_sub_path = Tuple{Vector{String}, Int, String}

        path  = ac_path.path
        len_p = ac_path.len_path
        summ  = d[path[1]]
        mtn_node = nothing

        for i in 2:len_p
            node = path[i]
            summ += d[node]
            if node in L_M_set && y_val[node] >= 0.9
                mtn_node = node
            end
            if summ > max_flt
                mtn_sub_path = (path[1:i], i, mtn_node) #Changer ceci pourait régler le problème
                break
            end  
        end
        return mtn_sub_path

    elseif result.status == "INFEASIBLE"
        infeas_sub_path = Tuple{Vector{String}, Int}
        
        path  = ac_path.path
        len_p = ac_path.len_path
        summ  = d[path[1]]

        for i in 2:len_p
            node = path[i]
            pred_node = path[i-1]
            summ += d[node]
            if (pred_node, node) in A_M_set
                break
            elseif summ > max_flt
                infeas_sub_path = (path[1:i], i)
                break
            end
        end
        return infeas_sub_path
    end
end



