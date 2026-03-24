
include("../structures.jl")

function verify_binarity(sx, sy, A, L_M)
    for a in A 
        if 0 < sx[a] && sx[a] < 1 
            println("IL S'AGIT DE X ")
            println(sx[a])
            return false
            break
        end 
    end  
    for j in L_M
        if 0 < sy[j] && sy[j] < 1 
            println("IL S'AGIT DE Y ")
            println(sy[j])
            return false
            break
        end 
    end 
    return true
end 

function buildM(env, instance_data, nbr_thread, silent, graph_reduc, time_limit)
    if graph_reduc
        instance_data = graph_reduction(instance_data)
    end

    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data
    other_data = instance_data.other_data

    nbr_K = instance_data.nbr_ac
    max_flt = instance_data.max_flying_time

    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes
    L_M = node_sets.mtn_nodes
    H_aircraft = node_sets.H_aircraft
    V = graph.nodes
    V_wt_st = vcat(a_nodes, fl_nodes)  #Nodes without st nodes

    A_S = arc_sets.arcs_S
    A_T = arc_sets.arcs_T
    A_M = arc_sets.arcs_M
    A_M_bar = arc_sets.arcs_M_bar
    A = graph.arcs

    f = fl_data.init_flying_time
    d = fl_data.d
   
    if graph_reduc
        F_bar = other_data["F_bar"]
        d_bar = other_data["d_bar"]
    end

    predecessors, successors = Dict(), Dict()
    predecessors, successors = update_neighborhood(A, predecessors, successors)
    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(A, pred_A_M, succ_A_M)
    
    # ===================== Model =====================
    model = Model(() -> Gurobi.Optimizer(env))
    #GRBsetintparam(env, "OutputFlag", 0)
    if silent
        set_optimizer_attribute(model, "OutputFlag", 0)
    end 
    set_optimizer_attribute(model, "Threads", nbr_thread)
    set_optimizer_attribute(model, "TimeLimit", time_limit) 
    set_optimizer_attribute(model, "MIPGap", 0.05) #5% de gap relatif
    # Gap absolu de 5
    #set_optimizer_attribute(model, "MIPGapAbs", 30.0)
    #set_optimizer_attribute(model, "LogFile", path_file)
    
    
    ft_obj = 0.0
    # ===================== Decision variables =====================
    model[:x] = @variable(model, 0 <= x[k in A] <= 1)                          #Arc (i,j) selected
    model[:y] = @variable(model,  0 <= y[j in L_M] <= 1)                        #1 if a maintenance is performed at node j 

    # ===================== Route Constraints =====================
    #= @constraint(model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(model, sum(x[i] for i in A_T) == nbr_K)
     =#
    @constraint(model, c1[i in V_wt_st], sum(x[(i,j)] for j in get(successors, i, [])) == 1)
    @constraint(model, c2[i in V_wt_st], sum(x[(j,i)] for j in get(predecessors, i, [])) == sum(x[(i,j)] for j in get(successors, i, [])))
    @constraint(model, c3[j in L_M], y[j] <= sum(x[(i, j)] for i in get(pred_A_M, j, [])))
    @constraint(model, sum(y[j] for j in L_M) <= length(H_aircraft))
     
    # =========== Variables ===========                        
    if graph_reduc 
        model[:u] = @variable(model, d_bar[j] <= u[j in V] <= F_bar[j])                  #Accumulatad flying time at node j 
    else
        model[:u] = @variable(model, d[j] <= u[j in V] <= max_flt)                  #Accumulatad flying time at node j 
    end 
    model[:rho] = @variable(model, rho[j in L_M] >= 0)                              #Remaining flying time at node j
    # =========== Objective ===========
    gamma_f = 1      
    ft_obj = sum(gamma_f*rho[j] for j in L_M)

    # ========  Flying time constraints ========
    @constraint(model, c4[(i,j) in A_M], rho[j] >= max_flt*x[(i, j)] - u[i] - (max_flt - d[i])*(1 - y[j]))
    @constraint(model, [(i,j) in A; j != "t"], u[j] <= u[i] + d[j] + (max_flt - d[i] - d[j])*(1 - x[(i,j)]), base_name = "c5[($i,$j)]")
    @constraint(model, c6[j in L_M], u[j] <= max_flt -(max_flt - d[j])*y[j])
    @constraint(model, [(i,j) in A_M_bar; j != "t"], u[j] >= u[i] + d[j] - max_flt*(1 - x[(i,j)]), base_name = "c7[($i,$j)]")
    @constraint(model, c8[(i,j) in A_M], u[j] >= u[i] + d[j] - max_flt*(1 - x[(i,j)]) - max_flt*y[j])
    @constraint(model, u["s"] == 0)
    @constraint(model, c9[k in a_nodes], u[k] == f[k])  

    @objective(model, Min, ft_obj)

    return model
end 

function solve_model(model, instance_data, fix_set, mip_set, relax_set, sx, sy, rf_or_fo)
    graph = instance_data.graph
    node_sets = graph.node_sets
    fl_data = instance_data.fl_data
    other_data = instance_data.other_data

    st_nodes = node_sets.dummy_nodes
    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes
    L_M = node_sets.mtn_nodes
    V = graph.nodes
    V_wt_st = vcat(a_nodes, fl_nodes)

    A = graph.arcs

    arc_day = other_data["arc_day"]
    
    x = model[:x]
    y = model[:y]
    u = model[:u]
    rho = model[:rho]
    
    # IMPORTANT: D'abord tout réinitialiser à l'état continu de base
    for arc in A
        # Unfixer si fixé
        if is_fixed(x[arc])
            unfix(x[arc])
        end
        # Enlever la contrainte binaire si elle existe
        if is_binary(x[arc])
            unset_binary(x[arc])
        end
    end
    for j in L_M
        # Unfixer si fixé
        if is_fixed(y[j])
            unfix(y[j])
        end
        # Enlever la contrainte binaire si elle existe
        if is_binary(y[j])
            unset_binary(y[j])
        end
    end 
    # Fixer les variables du fix_set
    L_M_set = Set(L_M)

    for day in fix_set
        arcs_set = [arc for (arc, v) in pairs(arc_day) if v == day]
        fl_set = Set(arc[2] for arc in arcs_set if arc[2] in L_M_set)
        
        for arc in arcs_set
            fix(x[arc], sx[arc]; force = true)
        end
        for j in fl_set
            fix(y[j], sy[j]; force = true)
        end
    end

    # Rendre binaires les variables du mip_set
    for day in mip_set    
        arcs_set = [arc for (arc, v) in pairs(arc_day) if v == day] 
        fl_set = Set(arc[2] for arc in arcs_set if arc[2] in L_M_set)
        
        for arc in arcs_set
            set_binary(x[arc])
        end
        for j in fl_set
            set_binary(y[j])
        end 
    end
    
    # Les variables du relax_set restent continues (état par défaut après unset_binary)
    # Pas besoin de code supplémentaire pour relax_set
    
    write_to_file(model, "modele.lp")
    optimize!(model)
    
    obj = objective_value(model)
    sx = Dict((i,j) => round(JuMP.value(x[(i,j)]), digits = 3) for (i,j) in A)
    sy = Dict(j => round(JuMP.value(y[j]), digits = 3) for j in L_M)
    su = Dict(j => round(JuMP.value(u[j]), digits = 3) for j in V)
    s_rho = Dict(j => round(Int, JuMP.value(rho[j])) for j in L_M)
    
    return (obj = obj, sx = sx, sy = sy, su = su, s_rho = s_rho, model = model)
end
 
function relax_and_fix(model, instance_data, output_file, size_day, overlap)
    verbose = false
    graph = instance_data.graph
    node_sets = graph.node_sets
    other_data = instance_data.other_data

    L_M = node_sets.mtn_nodes
    L_MS = node_sets.mtn_nodes_stn 
    MS = other_data["mtn_stations"]
    V = graph.nodes
    A = graph.arcs
    
    arc_day = other_data["arc_day"]
    all_days = sort(unique(values(arc_day)))
    sx = Dict{Tuple{String,String}, Float64}()
    sy = Dict{String, Float64}()

    begin_time = time()
    start_day = minimum(all_days)
    end_day = maximum(all_days)
    fix_set = Int[]
    mip_set = [i for i in start_day:start_day + size_day-1]
    relax_set = [i for i in start_day + size_day:end_day]
    rf_or_fo = "RF"
    iter = 1
    
    result = solve_model(model, instance_data, fix_set, mip_set, relax_set, sx, sy, rf_or_fo)
    curseur = start_day
    no_overlap = size_day - overlap
    sx = result.sx
    sy = result.sy
    model = result.model
    mdl = result.model

    write_both(output_file, "\tItération : $iter")
    write_both(output_file, "\tobjective_value = $(result.obj)\n") 

    if verbose
        println("\n\n--- ITERATION ", iter, " ---")
        println("\nfix_set = ", fix_set)
        println("\nmip_set = ", mip_set)
        println("\nrelax_set = ", relax_set)  
        println("\nno_overlap = ", no_overlap)  
    end

    while curseur + size_day-1 < end_day 
        iter += 1
        #mise à jour des ensembles
        curseur += no_overlap
        #println("\nCurrent curseur at day ", curseur)
        mip_set = [i for i in curseur:min(curseur + size_day-1, end_day)]
        relax_set = [i for i in curseur+size_day:end_day]
        fix_set = [i for i in start_day:curseur - 1]

        if verbose 
            println("\n\n--- ITERATION ", iter, " ---")
            println("\nfix_set = ", fix_set)
            println("\nmip_set = ", mip_set)
            println("\nrelax_set = ", relax_set)  
        end 
        result = solve_model(model, instance_data, fix_set, mip_set, relax_set, sx, sy, rf_or_fo)
        sx = result.sx
        sy = result.sy
        obj = result.obj
        model = result.model
        write_both(output_file, "\tItération : $iter")
        write_both(output_file, "\tobjective_value = $obj\n")
        #if all(isinteger, sx)
        if all(isa(v, Integer) for v in values(sx)) 
            break
        end
    end
    timeElapsed = round(time() - begin_time, digits = 4) 
    sv, sw = Dict(j => 0 for j in V), Dict(j => 0 for j in V)
    s_lambda, s_phi = Dict(j => 0 for j in L_M), Dict(j => 0 for j in L_M)

    mtn_stations_used = []
    active_j = Set(j for j in L_M if result.sy[j] >= 0.9)
    remaining_stations = setdiff(Set(MS), mtn_stations_used)

    for ms in remaining_stations
        if !isempty(intersect(Set(L_MS[ms]), active_j))
            push!(mtn_stations_used, ms)
        end
    end
    
    nbr_sts_used = round(Int, length(mtn_stations_used))

    
    other_info = Dict("model" => model, "time" => timeElapsed, "status" => termination_status(model), 
                    "mtn_stations_used" => mtn_stations_used, "nbr_sts_used" => nbr_sts_used)
    solution = Solution(result.obj, result.sx, result.sy, result.su, sv, sw, result.s_rho, s_lambda, s_phi, other_info)
    #return (obj = result.obj, x = result.sx, y = result.sy, u = result.su, rho = result.s_rho, time = timeElapsed, iter = iter)
    
    write_both(output_file,"\tON VERIFIE SI X ET Y SONT BINAIRES : $(verify_binarity(sx, sy, A, L_M))")
    
    return solution
end 

function fix_and_optimize(model, sx, sy, instance_data, output_file, size_day, overlap)
    verbose = false
    graph = instance_data.graph
    node_sets = graph.node_sets
    other_data = instance_data.other_data

    L_M = node_sets.mtn_nodes
    L_MS = node_sets.mtn_nodes_stn 
    MS = other_data["mtn_stations"]
    V = graph.nodes
    A = graph.arcs
    
    arc_day = other_data["arc_day"]
    all_days = sort(unique(values(arc_day)))
    
    begin_time = time()
    start_day = minimum(all_days)
    end_day = maximum(all_days)
    mip_set = [i for i in start_day:start_day + size_day-1]
    fix_set = [i for i in start_day + size_day:end_day]
    curseur = start_day
    no_overlap = size_day - overlap
    result = NamedTuple()

    rf_or_fo = "RF"
    iter = 0
   
    while curseur + size_day-1 <= end_day 
        iter += 1
        relax_set = []
        if verbose
            println("\n\n--- ITERATION ", iter, " ---")
            println("\nfix_set = ", fix_set)
            println("\nmip_set = ", mip_set)
            println("\nrelax_set = ", relax_set)
        end 
        result = solve_model(model, instance_data, fix_set, mip_set, relax_set, sx, sy, rf_or_fo)
        sx = result.sx
        sy = result.sy
        obj = result.obj
        model = result.model
        write_both(output_file, "\tItération : $iter")
        write_both(output_file, "\tobjective_value = $obj\n")
        #if all(isinteger, sx)
        if all(isa(v, Integer) for v in values(sx)) 
            break
        end
        #mise à jour des ensembles
        curseur += no_overlap
        #println("Curseur = ", curseur)
        #println("\nCurrent curseur at day ", curseur)
        mip_set = [i for i in curseur:min(curseur + size_day-1, end_day)]
        fix_set = vcat([i for i in start_day:curseur - 1], [i for i in curseur+size_day:end_day]) 
    end
    timeElapsed = round(time() - begin_time, digits = 4) 
    sv, sw = Dict(j => 0 for j in V), Dict(j => 0 for j in V)
    s_lambda, s_phi = Dict(j => 0 for j in L_M), Dict(j => 0 for j in L_M)

    mtn_stations_used = []
    active_j = Set(j for j in L_M if result.sy[j] >= 0.9)
    remaining_stations = setdiff(Set(MS), mtn_stations_used)

    for ms in remaining_stations
        if !isempty(intersect(Set(L_MS[ms]), active_j))
            push!(mtn_stations_used, ms)
        end
    end
    
    nbr_sts_used = round(Int, length(mtn_stations_used))

    
    other_info = Dict("time" => timeElapsed, "status" => termination_status(model), 
                    "mtn_stations_used" => mtn_stations_used, "nbr_sts_used" => nbr_sts_used)
    solution = Solution(result.obj, result.sx, result.sy, result.su, sv, sw, result.s_rho, s_lambda, s_phi, other_info)
    #return (obj = result.obj, x = result.sx, y = result.sy, u = result.su, rho = result.s_rho, time = timeElapsed, iter = iter)
    
    write_both(output_file,"\tON VERIFIE SI X ET Y SONT BINAIRES : $(verify_binarity(sx, sy, A, L_M))")

    return solution

end 