
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

function build_aircraft_path(start_aircraft::String, succ::Dict{String, String}, su)
    """Construit le chemin complet d'un avion à partir du graphe des successeurs"""
    chemin = ["s", start_aircraft]
    u_value = [su["s"], su[start_aircraft]]
    current = start_aircraft
    
    while haskey(succ, current) && succ[current] != "t"
        current = succ[current]
        push!(chemin, current)
        push!(u_value, su[current])
    end
    if haskey(succ, current)
        push!(u_value, su[succ[current]])

    end
    # Ajouter le nœud terminal si accessible
    haskey(succ, current) && push!(chemin, succ[current])
    
    return chemin, round.(u_value, digits = 2)
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


function print_solution(solution, instance_data, output_file, heuristic, silent=false)
    other_info = solution.other_info
    status = other_info["status"]
    graph = instance_data.graph
    node_sets = graph.node_sets
    L_M = node_sets.mtn_nodes
    A = graph.arcs
    a_nodes = node_sets.ac_nodes

    # Vérifier si la solution est valide
    valid_statuses = (MOI.OPTIMAL, MOI.FEASIBLE_POINT, MOI.TIME_LIMIT, MOI.INTERRUPTED)
    !(status in valid_statuses) && return Dict{String, String}()
    
    # Extraire les données de l'instance et de la solution
    obj_val, time_used, sx, sy = solution.obj, other_info["time"], solution.x, solution.y
    rho = solution.rho
    su = solution.u

    # Afficher les résultats généraux
    if heuristic
        write_both(output_file, "✅ Solution trouvée")
    else
        write_both(output_file, "✅ Solution trouvée (status: $status)")
    end
    write_both(output_file,"🎯 Objective global: $obj_val")
    write_both(output_file,"⏱️ Time used: $time_used")
    write_both(output_file,"")
    
    write_both(output_file,"✈️ Paths for each aircraft:")
    # Extraire les chemins des avions et trouver le plus long
    aircraft_paths = Dict{String, Vector{String}}()
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
        chemin, u_val = build_aircraft_path(aircraft, succ, su)
        println()
        #println(chemin)
        aircraft_paths[aircraft] = chemin
        # Vérifier si c'est le chemin le plus long
        
        path_with_values = ["$(chemin[i])_(u=$(u_val[i]))" for i in 1:length(chemin)]
        path_str = join(path_with_values, " ➡️  ")
        write_both(output_file, "🛩️ Aircraft $aircraft path ($(length(chemin)-3) nodes): $path_str")
        write_both(output_file,"")
    end
    maintenance_flights = [j for j in L_M if sy[j]>= 0.9]
    maintenance_at_node = Dict{String, String}()
    maintenance_info = Dict{String, Float64}()

    for j in maintenance_flights
        aircraft = find_aircraft_for_maintenance(j, aircraft_paths)
        maintenance_at_node[j] = aircraft       
        rho_val = rho[j]
        maintenance_info[aircraft] = rho_val
        write_both(output_file,"\n🛠️ Maintenance at node $j by aircraft $aircraft | ρ=$rho_val;")
    end
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