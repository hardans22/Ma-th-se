include("functions_for_benders.jl")

function build_sp(env, instance_data)
    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data
    other_data = instance_data.other_data

    max_flt = instance_data.max_flying_time
    
    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes
    H_aircraft = node_sets.H_aircraft
    L_M = node_sets.mtn_nodes
    V = graph.nodes
    V_wt_st = vcat(a_nodes, fl_nodes) 

    A = graph.arcs
   
    A_M = arc_sets.arcs_M
    A_M_bar = arc_sets.arcs_M_bar

    f = fl_data.init_flying_time
    d = fl_data.d
    pred_A_M = other_data["pred_A_M"]
    mtn_stations = other_data["mtn_stations"]
    

    sp_model = Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(sp_model, "OutputFlag", 0)
    set_optimizer_attribute(sp_model, "InfUnbdInfo", 1)
    set_optimizer_attribute(sp_model, "DualReductions", 0)

    # ===================Decision variables ===================
    sp_model[:x_copy] = @variable(sp_model, x_copy[k in A], Bin)                          #Copy of x
    sp_model[:y] = @variable(sp_model, y[j in L_M], Bin)        
    sp_model[:u] = @variable(sp_model, d[j] <= u[j in V] <= max_flt)                #Accumulatad flying time at node j  
    sp_model[:rho] = @variable(sp_model, rho[j in L_M] >= 0)                            #Remaining flying time at node j

    # =========== Objective ===========
    @objective(sp_model, Min,  sum(rho[j] for j in L_M))

    # ========  Flying time constraints ========

    @constraint(sp_model, sum(y[j] for j in L_M) <= length(H_aircraft))
    @constraint(sp_model, c3[j in L_M], y[j] <= sum(x_copy[(i,j)] for i in get(pred_A_M, j, [])))
    
    @constraint(sp_model, c4[(i,j) in A_M], rho[j] >= max_flt*x_copy[(i, j)] - u[i] - (max_flt - d[i])*(1 - y[j]))
    @constraint(sp_model, c5[(i,j) in A; j != "t"], u[j] <= u[i] + d[j] + (max_flt - d[i] - d[j])*(1 - x_copy[(i,j)]))
    @constraint(sp_model, c6[j in L_M], u[j] <= max_flt -(max_flt - d[j])*y[j])
    @constraint(sp_model, c7[(i,j) in A_M_bar; j != "t"], u[j] >= u[i] + d[j] - max_flt*(1 - x_copy[(i,j)]))
    @constraint(sp_model, c8[(i,j) in A_M], u[j] >= u[i] + d[j] - max_flt*(1 - x_copy[(i,j)]) - max_flt*y[j])
    @constraint(sp_model, u["s"] == 0)
    @constraint(sp_model, c9[k in a_nodes], u[k] == f[k])

    return sp_model
end


function solve_sp(sp_model, x_val, L_M)
    # Fixer x_copy
    x_copy = sp_model[:x_copy]
    for k in keys(x_val)
        fix(x_copy[k], x_val[k])
    end

    optimize!(sp_model)
    sp_status = termination_status(sp_model)
    if sp_status == MOI.OPTIMAL
        u_val = value.(sp_model[:u])
        y_val = value.(sp_model[:y])
        rho_val = value.(sp_model[:rho])
        return (status = sp_status, obj = objective_value(sp_model), time = round(solve_time(sp_model), digits = 4), y = y_val, u = u_val, rho = rho_val)
    else 
        #= if termination_status(sp_model) == MOI.INFEASIBLE
            println("MODEL INFEASIBLE")
        end  =#
        return (status = sp_status, obj = Inf, y = Dict(j => 0.0 for j in L_M), time = round(solve_time(sp_model), digits = 4))
    end
end 


#A REVOIR 
function solve_sp_analytic(x_val, instance_data)
    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data
    other_data = instance_data.other_data

    max_flt = instance_data.max_flying_time
    
    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes
    H_aircraft = node_sets.H_aircraft
    L_M = node_sets.mtn_nodes
    V = graph.nodes
    V_wt_st = vcat(a_nodes, fl_nodes) 

    A = graph.arcs
   
    A_M = arc_sets.arcs_M
    A_M_bar = arc_sets.arcs_M_bar

    f = fl_data.init_flying_time
    d = fl_data.d
    pred_A_M = other_data["pred_A_M"]
    mtn_stations = other_data["mtn_stations"]
    status = "OPTIMAL"
    time = @elapsed begin
        succ = Dict{String, String}()
        for (i, j) in Set(A)
            x_val[(i, j)] >= 0.9 && (succ[i] = j)
        end
        u = Dict(i => 0 for i in V)
        y = Dict(j => 0 for j in L_M)
        rho = Dict(j => 0 for j in L_M)
        u["s"] = 0
        u["t"] = 0
        aircraft_paths = Dict{String, Vector{String}}()
        mtn_subpath = Dict{String, Vector{String}}()
        infeas_subpath = Dict{String, Vector{String}}()
        found = false
        for aircraft in a_nodes
            chemin = ["s", aircraft]
            u[aircraft] = f[aircraft]
            current = aircraft
            current_mtn = nothing
            pred_current_mtn = nothing
            while haskey(succ, current) && succ[current] != "t"
                node_succ = succ[current]
                if (current, node_succ) in A_M
                    current_mtn = node_succ
                    pred_current_mtn = current
                end
                #println(current, " = ", u[current])
                if u[current] + d[node_succ] > max_flt
                    #println(current_mtn)
                    if current_mtn != nothing
                        y[current_mtn] = 1 
                        rho[current_mtn] = max_flt - u[pred_current_mtn]
                        u[current] = d[node_succ]
                        idx = findfirst(==(current_mtn), chemin)
                        if idx != nothing
                            mtn_subpath[aircraft] = chemin[1:findfirst(==(current_mtn), chemin)]
                        else 
                            mtn_subpath[aircraft] = vcat(chemin, current_mtn)
                        end 
                    else
                        status = "INFEASIBLE"
                        infeas_subpath[aircraft] = vcat(chemin, node_succ)
                    end     
                else
                    u[node_succ] = u[current] + d[node_succ]
                end 
                current = node_succ
                push!(chemin, current)
            end
            aircraft_paths[aircraft] = chemin
        end
        obj = sum(rho[j] for j in L_M)
    end    
    return (mtn_subpaths = mtn_subpath, infeas_subpaths = infeas_subpath, paths = aircraft_paths, u = u, rho = rho, y = y, obj = obj, time = time)

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
    
        summ = 0
        if result.status == MOI.OPTIMAL
            for i in 1:len_p
                if path[i] in L_M && y_val[path[i]] >= 0.9
                    mtn_sub_path[aircraft] = (path[1:i], i)
                    break
                end 
            end
        elseif result.status == MOI.INFEASIBLE
            mtn_poss = false
            for i in 1:len_p
                if (path[i], path[i+1]) in A_M
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


