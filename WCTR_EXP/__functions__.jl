# Fonctions helper optimisées (inchangées)
function construct_routes(xlsx_file)
    # Trier par avion et heure de départ
    df_flights = DataFrame(XLSX.readtable(xlsx_file, "Data"))
    df_sorted = sort(df_flights, [:TAIL_NUMBER, :DEPARTURE_TIME])

    # Ajouter l'encodage des vols
    transform!(df_sorted, [:ORIGIN_AIRPORT, :DESTINATION_AIRPORT, :DEPARTURE_TIME, :ARRIVAL_TIME] => 
            ((o, d, dep, arr) -> string.(o, "_", d, "_", dep, "_", arr)) => :FLIGHT_CODE)

    # Calculer le cumul des heures de vol par avion
    transform!(groupby(df_sorted, :TAIL_NUMBER), 
        :AIR_TIME => cumsum => :CUMULATIVE_FLYING_TIME)

    # Construire les itinéraires avec l'évolution du cumul
    itineraries = combine(groupby(df_sorted, :TAIL_NUMBER)) do group
        # Créer les codes de vol avec le cumul
        flight_evolution = ["$(row.FLIGHT_CODE) ($(row.CUMULATIVE_FLYING_TIME))" 
                            for row in eachrow(group)]
        
        DataFrame(
            NBR_FLIGHTS = nrow(group),
            ITINERARY = join(flight_evolution, " → "),
            TOTAL_FLYING_TIME = maximum(group.CUMULATIVE_FLYING_TIME)
        )
    end
    return itineraries
end 

function write_both(output_file, msg)
    println(msg)
    write(output_file, "\n$msg")
end

function build_aircraft_path(start_aircraft::String, succ::Dict{String, String}, su, sv, sw)
    """Construit le chemin complet d'un avion à partir du graphe des successeurs"""
    chemin = ["s", start_aircraft]
    u_value = [su["s"], su[start_aircraft]]
    v_value = [sv["s"], sv[start_aircraft]]
    w_value = [sw["s"], sw[start_aircraft]]
    current = start_aircraft
    
    while haskey(succ, current) && succ[current] != "t"
        current = succ[current]
        push!(chemin, current)
        push!(u_value, su[current])
        push!(v_value, sv[current])
        push!(w_value, sw[current])
    end
    if haskey(succ, current)
        push!(u_value, su[succ[current]])
        push!(v_value, sv[succ[current]])
        push!(w_value, sw[succ[current]])
    end
    # Ajouter le nœud terminal si accessible
    haskey(succ, current) && push!(chemin, succ[current])
    
    return chemin, round.(u_value, digits = 2), round.(v_value, digits = 2), round.(w_value, digits = 2)
end 

function find_aircraft_for_maintenance(maintenance_node::String, aircraft_paths::Dict{String, Vector{String}})
    """Trouve l'avion qui effectue la maintenance au nœud donné"""
    for (aircraft, path) in aircraft_paths
        maintenance_node ∈ path && return aircraft
    end
    return "UNKNOWN"
end

function calculate_path_distance(path::Vector{String}, d::Dict)
    """Calcule la distance totale d'un chemin"""
    total_distance = 0.0
    for i in 3:(length(path)-2)
        current_node = path[i]
        next_node = path[i+1]
        
        # Vérifier si la distance existe dans le dictionnaire
        if haskey(d, next_node)
            total_distance += d[next_node]
        end
    end
    
    return total_distance
end


function print_solution(solution, instance_data, output_file, silent=false)
    other_info = solution.other_info
    status = other_info["status"]
    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    L_M = node_sets.mtn_nodes
    A = graph.arcs
    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes

    # Vérifier si la solution est valide
    valid_statuses = (MOI.OPTIMAL, MOI.FEASIBLE_POINT, MOI.TIME_LIMIT, MOI.INTERRUPTED)
    !(status in valid_statuses) && return Dict{String, String}()
    
    # Extraire les données de l'instance et de la solution
    obj_val, time_used, sx, sy = solution.obj, other_info["time"], solution.x, solution.y
    rho, lambda, phi = solution.rho, solution.lambda, solution.phi
    su, sv, sw = solution.u, solution.v, solution.w
    mtn_stations_used, nbr_sts_used = other_info["mtn_stations_used"], other_info["nbr_sts_used"]
    #obj_rho, obj_lambda, obj_phi = solution["obj_rho"], solution["obj_lambda"], solution["obj_phi"]

   #=  # Calculer les statistiques des vols
    fl_satisfied = setdiff(keys(succ), union(["s"], a_nodes))
    not_fl_satisfied = setdiff(fl_nodes, keys(succ))
    nbr_mtn = sum(sy) =#
    
    # Fonction helper pour l'écriture
    
    # Afficher les résultats généraux
    write_both(output_file, "✅ Solution trouvée (status: $status)")
    write_both(output_file,"🎯 Objective global: $obj_val")
    #= write_both(output_file,"🎯 Objective rho: $obj_rho")
    write_both(output_file,"🎯 Objective lambda: $obj_lambda")
    write_both(output_file,"🎯 Objective phi: $obj_phi")
     =#write_both(output_file,"⏱️ Time used: $time_used")
   #=  write_both(output_file,"Nombre de vols satisfaits: $(length(fl_satisfied))")
    write_both(output_file,"Nombre de vols non satisfaits: $(length(not_fl_satisfied))") =#
    write_both(output_file,"")
    
    write_both(output_file,"✈️ Paths for each aircraft:")
    # Extraire les chemins des avions et trouver le plus long
    aircraft_paths = Dict{String, Vector{String}}()
    longest_aircraft = ""
    longest_path = String[]
    longest_length = 0
    aircraft_mtn = Dict{String, String}()
    # Construire le graphe des successeurs en une seule passe
    succ = Dict{String, String}()
    for (i, j) in A
        if sx[(i, j)] >= 0.9
            succ[i] = j
            #println("  Arc selected: $(i) -> $(j)")
        end
    end
    for aircraft in a_nodes
        #println(succ)
        chemin, u_val, v_val, w_val = build_aircraft_path(aircraft, succ, su, sv, sw)
        println()
        #println(chemin)
        aircraft_paths[aircraft] = chemin
        # Vérifier si c'est le chemin le plus long
        
        path_with_values = ["$(chemin[i])_(u=$(u_val[i]), v=$(v_val[i]), w=$(w_val[i]))" for i in 1:length(chemin)]
        path_str = join(path_with_values, " ➡️  ")
        write_both(output_file, "🛩️ Aircraft $aircraft path ($(length(chemin)-3) nodes): $path_str")
        write_both(output_file,"")
    end
    maintenance_flights = [j for j in L_M if sum(sy[j] for k in a_nodes) >= 0.9]
    maintenance_at_node = Dict{String, String}()
    maintenance_info = Dict{String, Tuple{Float64,Float64,Float64}}()

    for j in maintenance_flights
        aircraft = find_aircraft_for_maintenance(j, aircraft_paths)
        maintenance_at_node[j] = aircraft       
        rho_val, lambda_val, phi_val = rho[j], lambda[j], phi[j]
        maintenance_info[aircraft] = (rho_val, lambda_val, phi_val)
        write_both(output_file,"\n🛠️ Maintenance at node $j by aircraft $aircraft | ρ=$rho_val; λ=$lambda_val; φ=$phi_val")
    end
end


function compute_indicators(solution, instance_data, FH, FC, DY)
    other_info = solution.other_info
    status = other_info["status"]
    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data

    L_M = node_sets.mtn_nodes
    V = graph.nodes
    A = graph.arcs
    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes

    f = fl_data.init_flying_time
    h = fl_data.init_takeoff
    g = fl_data.init_flying_day
    b = fl_data.b
    d = fl_data.d
    tk = fl_data.tk

    max_flt = instance_data.max_flying_time
    max_tk = instance_data.max_takeoff
    max_day = instance_data.max_flying_day

    su, sv, sw, sx, sy = Dict(j => 0 for j in V), Dict(j => 0 for j in V), Dict(j => 0 for j in V), solution.x, solution.y
    s_rho, s_lambda, s_phi = Dict(j => 0 for j in L_M), Dict(j => 0 for j in L_M), Dict(j => 0 for j in L_M)
    
    succ = Dict{String, String}()
    for (i, j) in A
        sx[(i, j)] >= 0.9 && (succ[i] = j)
    end
    aircraft_starts = [j for (i, j) in A if i == "s" && sx[(i, j)] >= 0.9]
    aircraft_paths = Dict{String, Vector{String}}()

    for avion in a_nodes
        chemin, u_val, v_val, w_val = build_aircraft_path(avion, succ, su, sv, sw)
        aircraft_paths[avion] = chemin
        su[avion], sv[avion], sw[avion] = f[avion], h[avion], g[avion]
    end
    
    for path in values(aircraft_paths)
        for j in path
            if j == "s" || j == "t" || succ[j] == "t"
                continue
            end
            if succ[j] in L_M && sy[succ[j]] >= 0.9
                s_rho[succ[j]] = max_flt - su[j]
                s_lambda[succ[j]] = max_tk - sv[j]
                s_phi[succ[j]] = max_day - sw[j]
                su[succ[j]] = d[succ[j]]
                sv[succ[j]] = tk[succ[j]]
                sw[succ[j]] = 1
            else
                s_rho[succ[j]] = 0
                s_lambda[succ[j]] = 0
                s_phi[succ[j]] = 0
                su[succ[j]] = su[j] + d[succ[j]]
                sv[succ[j]] = sv[j] + tk[succ[j]]
                sw[succ[j]] = sw[j] + b[(j, succ[j])]
            end
        end
    end 

    feasible = 1
    for i in L_M
        if s_rho[i] < 0 || s_lambda[i] < 0 || s_phi[i] < 0
            feasible = 0
            break
        end   
    end
    solution.other_info["feasible"] = feasible
    active_j = Set(j for j in L_M if sy[j] >= 0.9)

    if FH
        solution.v, solution.w, solution.lambda, solution.phi = sv, sw, s_lambda, s_phi 
        obj_lambda = isempty(active_j) ? 0 : sum(s_lambda[j] for j in active_j)
        obj_phi = isempty(active_j) ? 0 : sum(s_phi[j] for j in active_j)
        solution.other_info["obj_lambda"], solution.other_info["obj_phi"] = obj_lambda, obj_phi
    end
    if FC
        solution.u, solution.w, solution.rho, solution.phi = su, sw, s_rho, s_phi
        obj_phi = isempty(active_j) ? 0 : sum(s_phi[j] for j in active_j)
        obj_rho = isempty(active_j) ? 0 : sum(s_rho[j] for j in active_j)
        solution.other_info["obj_rho"], solution.other_info["obj_phi"] = obj_rho, obj_phi
    end
    if DY
        solution.u, solution.v, solution.rho, solution.lambda = su, sv, s_rho, s_lambda
        obj_rho = isempty(active_j) ? 0 : sum(s_rho[j] for j in active_j)
        obj_lambda = isempty(active_j) ? 0 : sum(s_lambda[j] for j in active_j)
        solution.other_info["obj_rho"], solution.other_info["obj_lambda"] = obj_rho, obj_lambda
    end

    return solution
end 

function verification(solution, instance)
    A = solution["instance"]["A"]
    a_nodes = solution["instance"]["a_nodes"]
    L_M = solution["instance"]["L_M"]
    x_val = solution["x"]
    y_val = solution["y"]
    for a in A
        if sum(x_val[a,k] for k in a_nodes) > 1.01
            println("\nUN ARC EST AFFECTÉ À PLUS D'UN AVION\n")
        end
    end 
    for j in L_M
        if sum(y_val[j,k] for k in a_nodes) > 1.01
            println("UN NOEUD DE MAINTENANCE EST AFFECTÉ À PLUS D'UN AVION\n")
        end
    end 
end 