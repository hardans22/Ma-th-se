function build_aircraft_path(start_aircraft::String, succ::Dict{String, String})
    """Construit le chemin complet d'un aircraft à partir du graphe des successeurs"""
    path    = Vector{String}()
    push!(path, "s")
    push!(path, start_aircraft)
    current = start_aircraft
    len_p   = 2
    while haskey(succ, current) && succ[current] != "t"
        current = succ[current]
        push!(path, current)
        len_p += 1
    end
    len_p += 1
    # Ajouter le nœud terminal si accessible
    haskey(succ, current) && push!(path, succ[current])
    return path, len_p
end

function build_all_paths(x_val::Dict{Tuple{String, String}, Float64}, A::Vector{Tuple{String, String}}, a_nodes::Vector{String})
    """Construit les chemins de tous les avions à partir des variables x"""
    aircraft_paths = Dict{String, AircraftPath}()

    for aircraft in a_nodes
        # ✅ Un dictionnaire succ par avion
        succ = Dict{String, String}()
        for (i, j) in Set(A)
            x_val[(i, j)] >= 0.9 && (succ[i] = j)
        end

        path, len_p = build_aircraft_path(aircraft, succ)
        aircraft_paths[aircraft] = AircraftPath(path, len_p)
    end

    return aircraft_paths
end

function find_aircraft_for_maintenance(maintenance_node::String, aircraft_paths::Dict{String, AircraftPath})
    """Trouve l'aircraft qui effectue la maintenance au nœud donné"""
    for (aircraft, aircraftpath) in aircraft_paths
        path = aircraftpath.path
        maintenance_node in path && return aircraft
    end
    return "UNKNOWN"
end

function print_x_y(x_val, y_val, rho_val, A, L_M, d, H_aircraft)

    println("\nAFFICHAGE DE LA SOLUTION ")
    # Construire le graphe des successeurs en une seule passe
    succ = Dict{String, String}()
    for (i, j) in A
        x_val[(i, j)] >= 0.9 && (succ[i] = j)
    end
    
    # Extraire les chemins des aircrafts et trouver le plus long
    aircraft_paths = Dict{String, AircraftPath}()
    longest_aircraft = ""
    longest_path = String[]
    longest_length = 0
    
    # Identifier les aircrafts et construire leurs chemins
    #aircraft_starts = [("s", j) for j in a_nodes if x_val[("s", j)] >= 0.9]
    
    println("✈️ Paths for each aircraft:")
    aircraft_mtn = Dict{String, Vector{String}}()
    for aircraft in H_aircraft
        chemin, len_sp = build_aircraft_path(aircraft, succ)
        aircraft_paths[aircraft] = AircraftPath(chemin, len_sp)
        
        cumul = 0
        path_with_cumul = String[]
        for i in 1:length(chemin)
            cumul += d[chemin[i]]
            push!(path_with_cumul, "$(chemin[i])_(d=$(d[chemin[i]]), Σ=$cumul)")
        end

        path_str = join(path_with_cumul, " ➡️  ")
        println("\n🛩️ Aircraft $aircraft path ($(length(chemin)-3) nodes): $path_str")
        println("   📊 Cumul total des d: $cumul")
        
    end
    maintenance_flights = [j for j in L_M if y_val[j] >= 0.9]
    maintenance_at_node = Dict{String, String}()
    maintenance_info = Dict{String, Float64}()

    for j in maintenance_flights
        aircraft = find_aircraft_for_maintenance(j, aircraft_paths)
        maintenance_at_node[j] = aircraft       
        rho_v = rho_val[j]
        maintenance_info[aircraft] = rho_v
        println("\n🛠️ Maintenance at node $j by aircraft $aircraft | ρ=$rho_v;")
    end
    
end


function print_u(u, aircraft_paths)
    for (aircraft, aircraftpath) in aircraft_paths
        path = aircraftpath.path
        println("Aircraft: $aircraft")
        println(join([string(round(u[x], digits=2)) for x in path], " → "))
        println()
    end 
end 



# ========================================================

#Les fonctions d'affichage utilisé dans le __main__

# ========================================================

function write_both(output_file, msg)
    println(msg)
    write(output_file, "\n$msg")
end

function build_path(start_aircraft::String, succ::Dict{String, String}, su)
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

function find_aircraft_maintenance(maintenance_node::String, aircraft_paths::Dict{String, Vector{String}})
    """Trouve l'avion qui effectue la maintenance au nœud donné"""
    for (aircraft, path) in aircraft_paths
        maintenance_node ∈ path && return aircraft
    end
    return "UNKNOWN"
end


function print_solution(solution, instance_data, output_file, silent=false)
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
        chemin, u_val = build_path(aircraft, succ, su)
        println()
        #println(chemin)
        aircraft_paths[aircraft] = chemin
        # Vérifier si c'est le chemin le plus long
        
        path_with_values = ["$(chemin[i])_(u=$(u_val[i]))" for i in 1:length(chemin)]
        path_str = join(path_with_values, " ➡️  ")
        write_both(output_file, "🛩️ Aircraft $aircraft path ($(length(chemin)-3) nodes): $path_str")
        write_both(output_file,"")
    end
    maintenance_flights = [j for j in L_M if sum(sy[j] for k in a_nodes) >= 0.9]
    maintenance_at_node = Dict{String, String}()
    maintenance_info = Dict{String, Float64}()

    for j in maintenance_flights
        aircraft = find_aircraft_maintenance(j, aircraft_paths)
        maintenance_at_node[j] = aircraft       
        rho_val = rho[j]
        maintenance_info[aircraft] = rho_val
        write_both(output_file,"\n🛠️ Maintenance at node $j by aircraft $aircraft | ρ=$rho_val;")
    end
end
