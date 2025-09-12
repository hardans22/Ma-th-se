function print_solution(solution, output_file, silent=false)
    instance, status = solution["instance"], solution["status"]
    
    # Vérifier si la solution est valide
    valid_statuses = (MOI.OPTIMAL, MOI.FEASIBLE_POINT, MOI.TIME_LIMIT, MOI.INTERRUPTED)
    status ∉ valid_statuses && return Dict{String, String}()
    
    # Extraire les données de l'instance et de la solution
    A, L_M, a_nodes, fl_nodes = instance["A"], instance["L_M"], instance["a_nodes"], instance["fl_nodes"]
    obj_val, time_used, sx, sy = solution["obj"], solution["time"], solution["x"], solution["y"]
    d = instance["d"]
    rho, lambda, phi = solution["rho"], solution["lambda"], solution["phi"]
    su, sv, sw = solution["u"], solution["v"], solution["w"]
    mtn_stations_used, nbr_sts_used = solution["mtn_stations_used"], solution["nbr_sts_used"]
    
    # Construire le graphe des successeurs en une seule passe
    succ = Dict{String, String}()
    for (i, j) in A
        sx[(i, j)] >= 0.9 && (succ[i] = j)
    end
    
    # Calculer les statistiques des vols
    fl_satisfied = setdiff(keys(succ), union(["s"], a_nodes))
    not_fl_satisfied = setdiff(fl_nodes, keys(succ))
    nbr_mtn = sum(sy)
    
    # Fonction helper pour l'écriture
    function write_both(msg)
        !silent && println(msg)
        write(output_file, "\n$msg")
    end
    
    # Afficher les résultats généraux
    write_both("✅ Solution trouvée (status: $status)")
    write_both("🎯 Objective value: $obj_val")
    write_both("⏱️ Time used: $time_used")
    write_both("Nombre de vols satisfaits: $(length(fl_satisfied))")
    write_both("Nombre de vols non satisfaits: $(length(not_fl_satisfied))")
    write_both("")
    
    # Extraire les chemins des avions et trouver le plus long
    aircraft_paths = Dict{String, Vector{String}}()
    longest_aircraft = ""
    longest_path = String[]
    longest_length = 0
    
    # Identifier les avions et construire leurs chemins
    aircraft_starts = [(i, j) for (i, j) in A if i == "s" && sx[(i, j)] >= 0.9]
    
    write_both("✈️ Paths for each aircraft:")
    
    for (_, avion) in aircraft_starts
        chemin, u_value, v_value, w_value = build_aircraft_path(avion, succ, su, sv, sw)
        aircraft_paths[avion] = chemin
        
        # Vérifier si c'est le chemin le plus long
        if length(chemin) > longest_length
            longest_length = calculate_path_distance(chemin, d)
            longest_path = chemin
            longest_aircraft = avion
        end
        
        path_str = join(chemin, " ➡️  ")
        write_both("\n🛩️ Aircraft $avion path ($(length(chemin)) nodes): $path_str")
    end
    
    #= # Afficher la route la plus longue
    write_both("")
    write_both("🏆 ROUTE LA PLUS LONGUE:")
    write_both("✈️ Aircraft: $longest_aircraft")
    write_both("📏 Durée: $longest_length")
    write_both("🛤️ Chemin: $(join(longest_path, " ➡️  "))") =#
    
    # Identifier les maintenances actives et les mapper aux avions
    write_both("")
    maintenance_flights = [j for j in L_M if sy[j] >= 0.9]
    maintenance_at_node = Dict{String, String}()
    
    for j in maintenance_flights
        aircraft = find_aircraft_for_maintenance(j, aircraft_paths)
        maintenance_at_node[j] = aircraft
        
        rho_val, lambda_val, phi_val = value(rho[j]), value(lambda[j]), value(phi[j])
        write_both("🛠️ Maintenance at node $j by aircraft $aircraft | ρ=$rho_val; λ=$lambda_val; φ=$phi_val")
    end
    
    # Afficher les statistiques finales
    write_both("")
    write_both("Nombre de maintenances au total: $nbr_mtn")
    write_both("Stations de maintenance utilisées: $mtn_stations_used")
    write_both("Nombre de stations de maintenance utilisées: $nbr_sts_used")
    write_both("")
    
    # Retourner des informations sur la route la plus longue
    return Dict(
        "maintenance_at_node" => maintenance_at_node,
        "longest_aircraft" => longest_aircraft,
        "longest_path" => longest_path,
        "longest_length" => longest_length
    )
end

# Fonctions helper optimisées (inchangées)
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
    
    # Ajouter le nœud terminal si accessible
    haskey(succ, current) && push!(chemin, succ[current])
    
    return chemin, u_value, v_value, w_value
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