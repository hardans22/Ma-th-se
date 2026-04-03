include("functions_for_benders.jl")

function build_sp(env, instance_data)
    graph       = instance_data.graph
    node_sets   = graph.node_sets
    arc_sets    = graph.arc_sets
    fl_data     = instance_data.fl_data
    other_data  = instance_data.other_data

    
    fl_nodes    = node_sets.fl_nodes
    a_nodes     = node_sets.ac_nodes
    H_aircraft  = node_sets.H_aircraft
    V_wt_st     = vcat(a_nodes, fl_nodes) 
    max_flt     = instance_data.max_flying_time
    A_M_bar     = arc_sets.arcs_M_bar

    L_M = node_sets.mtn_nodes
    V   = graph.nodes
    A   = graph.arcs
    A_M = arc_sets.arcs_M
    f   = fl_data.init_flying_time
    d   = fl_data.d

    pred_A_M = other_data["pred_A_M"]
    mtn_stations = other_data["mtn_stations"]
    

    sp_model = Model(() -> Gurobi.Optimizer(env))
    set_optimizer_attribute(sp_model, "OutputFlag", 0)
    set_optimizer_attribute(sp_model, "InfUnbdInfo", 1)
    set_optimizer_attribute(sp_model, "DualReductions", 0)

    # ===================Decision variables ===================
    sp_model[:x_copy] = @variable(sp_model, x_copy[k in A])                          #Copy of x
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
    for arc in keys(x_val)
        fix(x_copy[arc], abs(x_val[arc]); force = true)
    end
    
    #write_to_file(sp_model, "sp.lp")

    optimize!(sp_model)
    sp_status = termination_status(sp_model)
    if sp_status == MOI.OPTIMAL
        u_val = value.(sp_model[:u])
        y_val = value.(sp_model[:y])
        rho_val = value.(sp_model[:rho])
        return (status = "OPTIMAL", obj = objective_value(sp_model), time = round(solve_time(sp_model), digits = 4), y = y_val, u = u_val, rho = rho_val)
    else 
        #= if termination_status(sp_model) == MOI.INFEASIBLE
            println("MODEL INFEASIBLE")
        end  =#
        return (status = "INFEASIBLE", obj = Inf, 
                y = Dict(j => 0.0 for j in L_M), 
                time = round(solve_time(sp_model), digits = 4))
    end
end 


#A REVOIR 
function solve_sp_pd(aircraft_paths, instance_data, a_nodes)
    graph       = instance_data.graph
    node_sets   = graph.node_sets
    arc_sets    = graph.arc_sets
    fl_data     = instance_data.fl_data

    L_M = node_sets.mtn_nodes
    V   = graph.nodes
    A_M = arc_sets.arcs_M
    f   = fl_data.init_flying_time
    d   = fl_data.d

    #a_nodes     = node_sets.ac_nodes
    max_flt     = instance_data.max_flying_time
    A_M_set      = Set(A_M)                  
    infeasible = false

    time = @elapsed begin
        u   = Dict{String, Float64}(j => 0.0 for j in V)
        y   = Dict{String, Int}(j => 0 for j in L_M)
        rho = Dict{String, Float64}(j => 0.0 for j in L_M)
        obj = 0.0
        status = "OPTIMAL"
        for aircraft in a_nodes
            path = aircraft_paths[aircraft].path
            len_path = aircraft_paths[aircraft].len_path
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
        end 
        u["t"] = 0.0
    end 
    if infeasible
        return (status = status, obj = Inf, time = round(time, digits = 6))
    end 
    return (status = status, obj = obj, y = y, u = u, 
            rho = rho, time = round(time, digits = 4))
end 


function solve_sp_analytic(x_val, instance_data) 
    graph        = instance_data.graph
    node_sets    = graph.node_sets
    arc_sets     = graph.arc_sets
    fl_data      = instance_data.fl_data
    other_data   = instance_data.other_data
    max_flt      = instance_data.max_flying_time
    a_nodes      = node_sets.ac_nodes

    L_M = node_sets.mtn_nodes
    V   = graph.nodes
    A_M = arc_sets.arcs_M
    f   = fl_data.init_flying_time
    d   = fl_data.d

    A_M_set      = Set(A_M)                  
    status       = "OPTIMAL"

    time = @elapsed begin
        
        succ = Dict{String, String}()
        for (i, j) in graph.arcs
            get(x_val, (i, j), 0.0) >= 0.9 && (succ[i] = j)
        end

        u   = Dict{String, Float64}(j => 0.0 for j in V)
        y   = Dict{String, Int}(j => 0 for j in L_M)
        rho = Dict{String, Float64}(j => 0.0 for j in L_M)
        aircraft_paths   = Dict{String, Vector{String}}()
        mtn_sub_path     = Dict{String, Tuple{Vector{String}, Int, String}}()
        infeas_sub_path  = Dict{String, Tuple{Vector{String}, Int}}()

        for aircraft in a_nodes
            chemin = String["s", aircraft]
            sizehint!(chemin, 50)          
            u[aircraft] = f[aircraft]
            current        = aircraft
            mtn_node       = nothing
            pred_mtn_node  = nothing
            idx            = 1

            while haskey(succ, current) && succ[current] != "t"
                node_succ = succ[current]

                if (current, node_succ) in A_M_set   
                    mtn_node      = node_succ
                    pred_mtn_node = current
                end

                if u[current] + d[node_succ] > max_flt
                    if mtn_node !== nothing           
                        y[mtn_node]   = 1
                        rho[mtn_node] = max_flt - u[pred_mtn_node]
                        u[current]    = d[node_succ]
                        mtn_sub_path[aircraft] = (chemin[1:idx], idx, mtn_node)
                        break
                    else
                        status = "INFEASIBLE"
                        path_copy = copy(chemin)      
                        push!(path_copy, node_succ)
                        infeas_sub_path[aircraft] = (path_copy, idx)
                        return(status = status, obj = Inf)
                    end
                else
                    u[node_succ] = u[current] + d[node_succ]
                end

                current = node_succ
                idx    += 1
                push!(chemin, current)
            end

            aircraft_paths[aircraft] = chemin
        end

        obj = sum(rho[j] for j in L_M)
    end

    return (mtn_subpaths = mtn_sub_path,
            infeas_subpaths = infeas_sub_path,
            paths = aircraft_paths, u = u,
            rho = rho, y = y, obj = obj,
            time = time, status = status)
end

function print_u(u, aircraft_paths)
    for (aircraft, aircraftpath) in aircraft_paths
        path = aircraftpath.path
        println("Aircraft: $aircraft")
        println(join([string(round(u[x], digits=2)) for x in path], " → "))
        println()
    end 
end 

function compare_solution(result, result_1, V, L_M)
    if result.status == MOI.OPTIMAL && result_1.status == "OPTIMAL"
        y_val   = result.y
        rho_val = result.rho
        u_val   = result.u

        y_val_1     = result_1.y
        rho_val_1   = result_1.rho
        u_val_1     = result_1.u
        eps = 1e-6
        for i in L_M
            if abs(u_val[i] - u_val_1[i]) > eps || abs(y_val[i] - y_val_1[i]) > eps || abs(rho_val[i] - rho_val_1[i]) > eps
                return false
            end
        end 
        for i in setdiff(V, L_M)
            if abs(u_val[i] - u_val_1[i]) > eps
                return false
            end
        end
    elseif result.status == result_1.status
        return true 
    else
        return false 
    end 
    return true
end 


#= function find_subpath_mtn(result, instance_data, aircraft_paths)
    graph      = instance_data.graph
    fl_data    = instance_data.fl_data
    A          = graph.arcs
    d          = fl_data.d
    max_flt    = instance_data.max_flying_time
    H_aircraft = graph.node_sets.H_aircraft

    L_M_set = Set(graph.node_sets.mtn_nodes)
    A_M_set = Set(graph.arc_sets.arcs_M)

    if result.status == "OPTIMAL"
        mtn_sub_path = Dict{String, Tuple{Vector{String}, Int}}()
        for aircraft in H_aircraft
            path  = aircraft_paths[aircraft]
            len_p = length(path)
            mtn_sub_path[aircraft] = (path[3:end], len_p-2) #Changer ceci pourait régler le problème
        end
        return mtn_sub_path

    elseif result.status == "INFEASIBLE"
        infeas_sub_path = Dict{String, Tuple{Vector{String}, Int}}()
        for aircraft in H_aircraft
            path  = aircraft_paths[aircraft]
            len_p = length(path)
            summ  = d[path[1]]

            for i in 2:len_p
                if (path[i-1], path[i]) in A_M_set
                    break
                elseif summ + d[path[i]] > max_flt
                    infeas_sub_path[aircraft] = (path[1:i], i)
                    break
                end
                summ += d[path[i]]
            end
        end
        return infeas_sub_path
    end

end 
 =#

function find_critical_subpath(instance_data, aircraft_paths)
    graph      = instance_data.graph
    fl_data    = instance_data.fl_data
    d          = fl_data.d
    max_flt    = instance_data.max_flying_time
    H_aircraft = graph.node_sets.H_aircraft

    ac_paths_reduced = Dict{String, AircraftPath}()

    for aircraft in H_aircraft
        path = aircraft_paths[aircraft].path
        len_p = aircraft_paths[aircraft].len_path
        sum = 0
        for i in 2:len_p
            sum += d[path[i]] 
            if sum > max_flt
                ac_paths_reduced[aircraft] = AircraftPath(path[1:i], i)
                break
            end 
        end 
    end 
    return ac_paths_reduced
end 

function find_best_subpath(result, instance_data, ac_paths_reduced)
    graph      = instance_data.graph
    fl_data    = instance_data.fl_data
    A          = graph.arcs
    d          = fl_data.d
    max_flt    = instance_data.max_flying_time
    H_aircraft = graph.node_sets.H_aircraft
    L_M_set    = Set(graph.node_sets.mtn_nodes)
    A_M_set    = Set(graph.arc_sets.arcs_M)

    if result.status == "INFEASIBLE"
        return ac_paths_reduced

    elseif result.status == "OPTIMAL"
        ref_obj = result.obj

        for aircraft in keys(ac_paths_reduced)
            current_path = copy(ac_paths_reduced[aircraft].path)

            while length(current_path) > 1
                candidate_path = current_path[1:end-1]

                # Test avec le chemin réduit
                ac_paths_reduced[aircraft].path = candidate_path
                rslt = solve_sp_pd(ac_paths_reduced, instance_data, keys(ac_paths_reduced))

                if abs(rslt.obj - ref_obj) >= 1e-6
                    # Ce retrait change l'objectif → current_path est le minimal utile
                    ac_paths_reduced[aircraft].path = current_path
                    break
                end

                # Pas de changement → on continue de réduire
                current_path = candidate_path
            end
        end
    end

    return ac_paths_reduced
end


function find_subpath_mtn(result, instance_data, aircraft_paths)
    graph      = instance_data.graph
    fl_data    = instance_data.fl_data
    A          = graph.arcs
    d          = fl_data.d
    max_flt    = instance_data.max_flying_time
    H_aircraft = graph.node_sets.H_aircraft

    L_M_set = Set(graph.node_sets.mtn_nodes)
    A_M_set = Set(graph.arc_sets.arcs_M)

    if result.status == "OPTIMAL"
        y_val = result.y
        rho_val = result.rho
        mtn_sub_path = Dict{String, Tuple{Vector{String}, Int, String}}()

        for aircraft in H_aircraft
            path  = aircraft_paths[aircraft].path
            len_p = aircraft_paths[aircraft].len_path
            summ  = d[path[1]]
            mtn_node = nothing

            for i in 2:len_p
                node = path[i]
                pred_node = path[i-1]
                if node in L_M_set && y_val[node] >= 0.9
                    #val = result.rho[node].
                    mtn_node = node
                end
                if summ + d[node] > max_flt
                    #println(mtn_node)
                    mtn_sub_path[aircraft] = (path[1:i], i, mtn_node) #Changer ceci pourait régler le problème
                    break
                    #= if mtn_node !== nothing && rho_val[mtn_node] == 0.0
                        mtn_sub_path[aircraft] = ([], 0, mtn_node) #Changer ceci pourait régler le problème
                        break   

                    else 
                        mtn_sub_path[aircraft] = (path[1:i], i, mtn_node) #Changer ceci pourait régler le problème
                        break
                    end =#
                end  
                summ += d[node]
            end
        end
        return mtn_sub_path

    elseif result.status == "INFEASIBLE"
        infeas_sub_path = Dict{String, Tuple{Vector{String}, Int}}()
        for aircraft in H_aircraft
            path  = aircraft_paths[aircraft].path
            len_p = aircraft_paths[aircraft].len_path
            summ  = d[path[1]]

            for i in 2:len_p
                if (path[i-1], path[i]) in A_M_set
                    break
                elseif summ + d[path[i]] > max_flt
                    infeas_sub_path[aircraft] = (path[1:i], i)
                    break
                end
                summ += d[path[i]]
            end
        end
        return infeas_sub_path
    end
end


function find_neighbor(path, len_p, A, d, max_flt)
    pred_A, succ_A = Dict(), Dict()
    pred_A, succ_A = update_neighborhood(A, pred_A, succ_A)
    list_path = Vector{Tuple{Vector{String}, Int}}()
    curr_node = path[end]
    sum_temp = sum(d[l] for l in path[1:end])
    for j in succ_A[curr_node]
        if sum_temp + d[j] > max_flt
            push!(list_path, (vcat(path,j), len_p+1))
        end 
    end 
    return list_path
end 


function verify_cuts(benders_cuts)
    for (key, cuts_list) in benders_cuts
        for cut in cuts_list
            lhs_val = value(cut.func)        # évalue l'expression symbolique
            rhs_val = cut.set.upper          # ou .lower selon le type de set
            slack = lhs_val - rhs_val
            
            if abs(slack) < 1e-6
                println("Cut $i actif")
            end
        end 
    end
end 


function check_cuts(benders_cuts, x_val; tol=1e-12)
    println("="^50)
    println("Vérification des Benders cuts")
    println("="^50)

    for (key, cuts_list) in benders_cuts
        for cut in cuts_list
            expr = cut.func  # AffExpr

            # Évaluer manuellement l'AffExpr
            lhs_val = expr.constant
            for (coef, var) in linear_terms(expr)
                if var in keys(x_val)  # ignore z et toute autre variable absente du dict
                    lhs_val += coef * x_val[var]
                end
            end

            rhs_val = cut.set.upper
            slack = rhs_val - lhs_val

            status = slack >= -tol ? "✅ satisfait" : "❌ violé"
            if slack < -0.0
                println("Coupe de l'iteration : $key : lhs = $(round(lhs_val, digits=4)) | rhs = $(round(rhs_val, digits=4)) | slack = $(round(slack, digits=4)) → $status")
                println(cut)
            end 
        end
    end 
    println("="^50)
end