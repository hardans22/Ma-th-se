using JSON
include("structures.jl")

# ====================================================================
# FONCTIONS UTILITAIRES POUR LA MANIPULATION DE GRAPHES
# ====================================================================

"""Retourne tous les arcs sortants du nœud i"""
delta_plus(i, A) = [(u, v) for (u, v) in A if u == i]

"""Retourne tous les arcs entrants vers le nœud i"""
delta_minus(i, A) = [(u, v) for (u, v) in A if v == i]

"""
Met à jour les ensembles d'arcs en remplaçant 'arc' par 'new_arc'
Utilise des Sets temporaires pour des recherches O(1)
"""
function arc_set_update(A, a_nodes, fl_nodes, flight_legs, init_airport, DT, AT, mtn_time, mtn_stations)
    A_K, A_F, A_M, L_M, L_a_M = [], [], [], [], []  # Arcs et nœuds de maintenance
    for x in Set(A)
        i, j = x
        # Détermination de la station d'arrivée
        if i in Set(fl_nodes)
            arrival_station = flight_legs[i][2]
            if j in Set(fl_nodes)
                push!(A_F, x)
            end
        elseif i in Set(a_nodes)
            arrival_station = init_airport[i]
            push!(A_K, x)
        end
        # Arc de maintenance possible ?
        if j in fl_nodes && arrival_station in mtn_stations && DT[j] - AT[i] >= mtn_time
            push!(A_M, x)
            push!(L_M, j)
            push!(L_a_M, i)
        end
    end
    unique!(L_M)
    unique!(L_a_M)
    #L_MS = Dict(ms => [i for i in L_M if flight_legs[i][1] == ms] for ms in mtn_stations)
    L_MS, M_FL_O, M_FL_D = Dict(ms => String[] for ms in mtn_stations), Dict(ms => String[] for ms in mtn_stations), Dict(ms => String[] for ms in mtn_stations)
    for i in L_M
        ms = flight_legs[i][1]
        haskey(L_MS, ms) && push!(L_MS[ms], i)
    end

    for i in fl_nodes
        ms_o, ms_d = flight_legs[i]
        haskey(M_FL_O, ms_o) && push!(M_FL_O[ms_o], i)
        haskey(M_FL_D, ms_d) && push!(M_FL_D[ms_d], i)
    end
    
    A_M_bar = setdiff(A, A_M)
    return A_K, A_F, A_M, A_M_bar, L_M, L_MS, M_FL_O, M_FL_D, L_a_M
end

"""
Met à jour les données (durées, temps) lors de la fusion de deux nœuds
"""
function data_update(d, tk, DT, AT, nodes_rm)
    # Suppression des anciens nœuds (sauf source et puits)
    for node in nodes_rm 
        for dict in [d, tk, DT, AT]
            pop!(dict, node, nothing)
        end
    end
    return d, tk, DT, AT
end

function fl_a_nodes_update(i, j, ij, a_nodes, fl_nodes, flight_legs, d, tk, init_airport, init_flying_time, init_takeoff, init_flying_day)
    flight_legs[ij] = []
    # Mise à jour des ensembles de nœuds
    if i in fl_nodes
        fl_nodes = setdiff(fl_nodes, [i])
        push!(flight_legs[ij], flight_legs[i][1])
        push!(flight_legs[ij], flight_legs[j][2])
        push!(fl_nodes, ij)
    elseif i in a_nodes
        a_nodes = setdiff(a_nodes, [i])
        push!(a_nodes, ij)
        init_airport[ij] = init_airport[i]
        # Transfert des propriétés d'avion
        init_flying_time[ij] = init_flying_time[i] + d[j]
        init_takeoff[ij] = init_takeoff[i] + tk[j]
        init_flying_day[ij] = init_flying_day[i]
        for dict in [init_flying_time, init_takeoff, init_flying_day]
            filter!(kv -> kv[1] != i, dict)
        end
    end
    fl_nodes = setdiff(fl_nodes, [j])
    
    return a_nodes, fl_nodes, flight_legs, init_airport, init_flying_time, init_takeoff, init_flying_day
end 

# Fonction pour mettre à jour les structures de voisinage
function update_neighborhood(A, predecessors, successors)
    empty!(predecessors)
    empty!(successors)
    for (u, v) in A
        if !haskey(predecessors, v)
            predecessors[v] = []
        end
        if !haskey(successors, u)
            successors[u] = []
        end
        push!(predecessors[v], u)
        push!(successors[u], v)
    end
    return predecessors, successors
end

# Fonction pour calculer b pour un arc spécifique
function calc_b(i, j, DT)
    return (floor(Int, DT[j]/1440) - floor(Int, DT[i]/1440)) 
end

# ====================================================================
# FONCTION PRINCIPALE DE CONSTRUCTION DU GRAPHE
# ====================================================================

"""
Construit le graphe de planification des vols avec preprocessing optionnel
- file: fichier JSON contenant l'instance
- preprocess: si true, applique la réduction de graphe
"""
function build_graph(file)
    instance = JSON.parsefile(file)
    # ================================================================
    # EXTRACTION DES DONNÉES D'ENTRÉE
    # ================================================================
    nbr_FL, nbr_K = instance["number_of_flight_legs"], instance["number_of_aircrafts"]
    f_legs, aircrafts = instance["flight_legs"], instance["aircrafts"]
    mtn_stations = Set(instance["maintenance_stations"])  # Set pour recherche O(1)
    init_airport = instance["initial_airport_aircraft"]
    TRT, mtn_time = instance["turn_around_time"], instance["maintenance_time"]
    init_flying_time = instance["initial_flying_time"]
    init_takeoff, init_flying_day = instance["initial_takeoff"], instance["initial_flying_day"]
    end_h_time, nbr_TP = instance["end_horizon_time"], instance["nbr_TP"]
    max_flt = instance["maximum_flying_time"]
    max_tk = instance["maximum_takeoff"]
    max_fd = instance["maximum_flying_day"]
    
    # ================================================================
    # CONSTRUCTION DES NŒUDS
    # ================================================================
    # Nœuds de vols : "origine_destination_départ_arrivée"
    fl_nodes = [f_legs[i][1]*"_"*f_legs[i][2]*"_"*string(f_legs[i][3])*"_"*string(f_legs[i][4]) 
                for i in 1:nbr_FL]
    a_nodes = aircrafts
    st_nodes = ["s", "t"]
    V = vcat(fl_nodes, a_nodes, st_nodes)
    #V_wt_st = vcat(a_nodes, fl_nodes)
    
    # ================================================================
    # INITIALISATION DES STRUCTURES DE DONNÉES
    # ================================================================
    # Dictionnaires de données des vols
    flight_legs = Dict(fl_nodes[i] => [f_legs[i][1], f_legs[i][2]] for i in 1:nbr_FL)
    DT = Dict(fl_nodes[i] => f_legs[i][3] for i in 1:nbr_FL)  # Departure Time
    AT = Dict(fl_nodes[i] => f_legs[i][4] for i in 1:nbr_FL)  # Arrival Time
    d = Dict(fl_nodes[i] => f_legs[i][5] for i in 1:nbr_FL)   # Duration
    tk = Dict(fl_nodes[i] => 1 for i in 1:nbr_FL)   # Takeoff
    
    H_aircraft = Vector{String}()
    # Initialisation des avions et nœuds source/puits
    temp = minimum(values(DT))
    A_S = [(st_nodes[1], k) for k in a_nodes]  # Arcs source -> avions
    for k in a_nodes
        d[k],tk[k] = init_flying_time[k], init_takeoff[k]
        if max_flt - init_flying_time[k] < 240
            DT[k], AT[k] = temp-mtn_time-35, temp-mtn_time-35  # Initialisation au plus tôt
        else
            DT[k], AT[k] = temp, temp  # Initialisation au plus tôt
        end
        if d[k] > 0
            push!(H_aircraft, k)
        end 
    end
    d["s"], tk["s"], DT["s"], AT["s"] = 0, 0, temp, temp 
    d["t"], tk["t"], DT["t"], AT["t"] = 0, 0, end_h_time, end_h_time
    
    # ================================================================
    # CONSTRUCTION DES ARCS
    # ================================================================
    A_K, A_F, A_T = [], [], []  # Arcs avions->vols, vols->vols, vols->puits
    #a = Dict()  # Matrice d'assignation temporelle
    a_day = Dict{String, Int}()
    println("Construction des arcs...")
    for i in fl_nodes
        origin, dest = flight_legs[i][1], flight_legs[i][2]
        dt_i, at_i = DT[i], AT[i]
        # A_K: connexions avion -> vol
        for k in a_nodes
            if init_airport[k] == origin && dt_i - AT[k] <= 1440
                push!(A_K, (k, i))
            end
        end
        # A_F: connexions vol -> vol (avec contrainte de temps de rotation)
        for j in Set(fl_nodes)
            if i != j && dest == flight_legs[j][1]
                time_diff = DT[j] - at_i
                if TRT <= time_diff <= 1440 #******************************************************************
                    push!(A_F, (i, j))
                end
            end
        end
        # A_T: connexions vol -> puits
        (end_h_time - at_i <= 1440) && push!(A_T, (i, st_nodes[2]))
        # Matrice d'assignation par période temporelle
        day = ceil(Int, dt_i / 1440)
        a_day[i] = day
    end
    
    # ================================================================
    # ARCS DE MAINTENANCE ET PARAMÈTRES BINAIRES
    # ================================================================
    A = vcat(A_S, A_K, A_F, A_T)
    A_M, L_M = [], []  # Arcs et nœuds de maintenance
    b = Dict()  # Variable binaire jour
    arc_info = Dict{Tuple{String,String},String}()
    println("Identification des opportunités de maintenance...")
    for x in Set(A)
        i, j = x
        # Détermination de la station d'arrivée
        arrival_station = if i in fl_nodes
            flight_legs[i][2]
        elseif i in a_nodes
            init_airport[i]
        end
        # Arc de maintenance possible ?
        if j in fl_nodes && arrival_station in mtn_stations && DT[j] - AT[i] >= mtn_time
            push!(A_M, x)
            push!(L_M, j)
            arc_info[x] = "A_M"
        else
            arc_info[x] = "A_M_bar"
        end
        # Variable binaire pour changement de jour
        b[x] = floor(Int, DT[j]/1440) - floor(Int, DT[i]/1440)
        if i == "s" && j in a_nodes
            b[x] = init_flying_day[j]  # Pas de changement de jour pour les arcs source->avion
        end
    end
    L_M = unique(L_M)
    A_M_bar = setdiff(A, A_M)
    
    L_MS, M_FL_O, M_FL_D = Dict(ms => String[] for ms in mtn_stations), Dict(ms => String[] for ms in mtn_stations), Dict(ms => String[] for ms in mtn_stations)

    for i in L_M
        ms = flight_legs[i][1]
        haskey(L_MS, ms) && push!(L_MS[ms], i)
    end

    for i in fl_nodes
        ms_o, ms_d = flight_legs[i]
        haskey(M_FL_O, ms_o) && push!(M_FL_O[ms_o], i)
        haskey(M_FL_D, ms_d) && push!(M_FL_D[ms_d], i)
    end

    
    NH_aircraft = setdiff(a_nodes, H_aircraft)
    #merge!(instance, Dict("M_FL_O" => M_FL_O, "M_FL_D" => M_FL_D, "a_day" => a_day))
    other_data = Dict("mtn_stations" => instance["maintenance_stations"], "fh_day" => instance["fh_day"], 
                    "fh_tk" => instance["fh_tk"], "M_FL_O" => M_FL_O, "M_FL_D" => M_FL_D, "a_day" => a_day, 
                    "ms_capacity" => instance["mtn_station_capacity"], "nbr_TP" => nbr_TP)

    node_sets = NodeSets(st_nodes, fl_nodes, a_nodes, L_M, L_MS, H_aircraft, NH_aircraft)
    arc_sets = ArcSets(A_S, A_K, A_F, A_T, A_M, A_M_bar)
    graph = Graph(V, A, node_sets, arc_sets, arc_info)
    flight_data = FlightData(flight_legs, d, b, tk, DT, AT, init_airport, init_flying_time, init_takeoff, init_flying_day)
    instance_data = InstanceData(graph, flight_data, other_data, TRT, mtn_time, max_flt, max_tk, max_fd, nbr_K)
    #println(L_MS)
    # Mise à jour de l'instance avec toutes les structures calculées
    
    #println("\nConstruction terminée. Réduction: $arc_reduc arcs et $node_reduc noeuds\n")
    return instance_data
end

function graph_reduction(instance_data::InstanceData)
    # ================================================================
    # PREPROCESSING OPTIMISÉ : RÉDUCTION DU GRAPHE
    # ================================================================
    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data
    other_data = instance_data.other_data

    mtn_time = instance_data.mtn_time
    
    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes
    
    A_S = arc_sets.arcs_S
    A_T = arc_sets.arcs_T
    A_M_bar = arc_sets.arcs_M_bar
    A = graph.arcs
    V = graph.nodes
    
    d = fl_data.d
    b = fl_data.b
    tk = fl_data.tk
    DT = fl_data.DT
    AT = fl_data.AT

    flight_legs = fl_data.fl_legs
    init_airport = fl_data.init_airport
    init_flying_time = fl_data.init_flying_time
    init_takeoff = fl_data.init_takeoff
    init_flying_day = fl_data.init_flying_day

    mtn_stations = other_data["mtn_stations"]

    old_arc_nbr = length(A)
    old_node_nbr = length(V)
    arc_reduc = 0
    node_reduc = 0
    println("\nNombre d'arcs avant : $(old_arc_nbr) arcs")
    println("Nombre de noeuds avant : $(old_node_nbr) noeuds")
    
    nodes_rm, nodes_new = [], []
    # Structures pour éviter les recalculs
    predecessors = Dict()
    successors = Dict()
    # Initialisation
    predecessors, successors = update_neighborhood(A, predecessors, successors)
    A_S_set = Set(A_S)
    A_T_set = Set(A_T)
    queue = copy(A_M_bar)
    visited = Set()
    while !isempty(queue)
        arc = popfirst!(queue)
        arc in visited && continue
        push!(visited, arc)
        i, j = arc
        merged = false
        # --------------------------------------------------------
        # CAS 1: i est l'unique prédécesseur de j
        # --------------------------------------------------------
        if !(arc in A_S_set) && get(predecessors, j, []) == [i]
            ij = i * "_" * j
            merged = true
            push!(nodes_new, ij)
            # Mise à jour des nœuds
            a_nodes, fl_nodes, flight_legs, init_airport, init_flying_time, init_takeoff, init_flying_day = 
                fl_a_nodes_update(i, j, ij, a_nodes, fl_nodes, flight_legs, d, tk, init_airport, init_flying_time, init_takeoff, init_flying_day)
            
            # Collecte des arcs à modifier
            inarc_i = [(u, i) for u in get(predecessors, i, [])]
            outarc_j = [(j, v) for v in get(successors, j, [])]
            outarc_i = [(i, v) for v in get(successors, i, [])]
            arcs_to_remove = vcat([(i, j)], inarc_i, outarc_j, outarc_i)
            # Nouveaux arcs
            # Mise à jour des données
            d[ij] = d[i] + d[j]
            tk[ij] = tk[i] + tk[j]
            DT[ij] = DT[i]
            AT[ij] = AT[j]
            # Mise à jour sélective de b : supprimer les anciens arcs et ajouter les nouveaux
            new_arcs = []
            for u in get(predecessors, i, [])
                new_arc = (u, ij)
                push!(new_arcs, new_arc)
                if DT[ij] - AT[u] <= mtn_time
                    push!(queue, new_arc)
                end
            end
            for v in get(successors, j, [])
                new_arc = (ij, v)
                push!(new_arcs, new_arc)
                if DT[v] - AT[ij] <= mtn_time
                    push!(queue, new_arc)
                end
            end
            
            for arc in Set(arcs_to_remove)
                delete!(b, arc)
            end
            for arc in Set(new_arcs)
                b[arc] = calc_b(arc[1], arc[2], DT)
            end
            # Mise à jour de A
            A = setdiff(A, Set(arcs_to_remove))
            append!(A, new_arcs)
            # Mise à jour des sets pour éviter les recalculs
            A_S_set = Set([(u, v) for (u, v) in A if u == "s"])
            A_T_set = Set([(u, v) for (u, v) in A if v == "t"])
            (i != "s") && push!(nodes_rm, i)
            (j != "t") && push!(nodes_rm, j)
            
        # --------------------------------------------------------
        # CAS 2: j est l'unique successeur de i
        # --------------------------------------------------------
        elseif !(arc in A_T_set) && get(successors, i, []) == [j]
            ij = i * "_" * j
            merged = true
            push!(nodes_new, ij)
            # Mise à jour des nœuds
            a_nodes, fl_nodes, flight_legs, init_airport, init_flying_time, init_takeoff, init_flying_day = 
                fl_a_nodes_update(i, j, ij, a_nodes, fl_nodes, flight_legs, d, tk, init_airport, init_flying_time, init_takeoff, init_flying_day)
            # Collecte des arcs à modifier
            inarc_i = [(u, i) for u in get(predecessors, i, [])]
            outarc_j = [(j, v) for v in get(successors, j, [])]
            inarc_j = [(u, j) for u in get(predecessors, j, [])]
            arcs_to_remove = vcat([(i, j)], inarc_i, outarc_j, inarc_j)
            # Nouveaux arcs
            # Mise à jour des données
            d[ij] = d[i] + d[j]
            tk[ij] = tk[i] + tk[j]
            
            DT[ij] = DT[i]
            AT[ij] = AT[j] 
            # Mise à jour sélective de b : supprimer les anciens arcs et ajouter les nouveaux
            new_arcs = []
            for u in get(predecessors, i, [])
                new_arc = (u, ij)
                push!(new_arcs, new_arc)
                if DT[ij] - AT[u] <= mtn_time
                    push!(queue, new_arc)
                end
            end
            for v in get(successors, j, [])
                new_arc = (ij, v)
                push!(new_arcs, new_arc)
                if DT[v] - AT[ij] <= mtn_time
                    push!(queue, new_arc)
                end
            end
            
            for arc in arcs_to_remove
                delete!(b, arc)
            end
            for arc in new_arcs
                b[arc] = calc_b(arc[1], arc[2], DT)
            end
            # Mise à jour de A
            A = setdiff(A, Set(arcs_to_remove))
            append!(A, new_arcs)
            # Mise à jour des sets
            A_S_set = Set([(u, v) for (u, v) in A if u == "s"])
            A_T_set = Set([(u, v) for (u, v) in A if v == "t"])
            (i != "s") && push!(nodes_rm, i)
            (j != "t") && push!(nodes_rm, j)
        end
        # Mise à jour du voisinage seulement si fusion effectuée
        if merged
            predecessors, successors = update_neighborhood(A, predecessors, successors)
            nodes_rm = []
            nodes_new = []
        end
    end
    
    # ============================================================
    # CALCULS FINAUX (UNE SEULE FOIS)
    # ============================================================
    # Mise à jour finale des structures de maintenance
    A_K, A_F, A_M, A_M_bar, L_M, L_MS, M_FL_O, M_FL_D, L_a_M = arc_set_update(A, a_nodes, fl_nodes, flight_legs, init_airport, DT, AT, mtn_time, mtn_stations)
    A_S, A_T = collect(A_S_set), collect(A_T_set)
    
    # Nettoyage final des données
    d, tk, DT, AT = data_update(d, tk, DT, AT, nodes_rm)
    V = vcat(["s", "t"], fl_nodes, a_nodes)
    for k in a_nodes
        d[k],tk[k] = init_flying_time[k], init_takeoff[k]
    end
    
    new_arc_nbr = length(A)
    new_node_nbr = length(V)
    arc_reduc = old_arc_nbr - new_arc_nbr
    node_reduc = old_node_nbr - new_node_nbr
    println("\nNombre d'arcs après : $new_arc_nbr arcs")
    println("Nombre de noeuds après : $new_node_nbr noeuds")
    # Recalcul de la matrice d'assignation temporelle
    #= a = Dict()
    for i in fl_nodes
        day = ceil(Int, DT[i] / 1440)
        for t in 1:nbr_TP
            a[(i, t)] = (day == t) ? 1 : 0
        end
    end =#
    # Initialiser le dictionnaire
    a_day = Dict{String, Int}()

    # Remplir le dictionnaire
    for i in fl_nodes
        day = ceil(Int, DT[i] / 1440)  # Adapter selon votre structure
        a_day[i] = day
    end
    other_data["a_day"] = a_day

    #= F_bar = Dict(j => instance_data.max_flying_time for j in fl_nodes)
    d_bar = Dict(j => d[j] for j in fl_nodes) =#
    
    # ============================================================
    # CALCULS DES BORNES : F_bar et d_bar
    # ============================================================
    predecessors = Dict()
    successors = Dict()
    # Initialisation
    predecessors, successors = update_neighborhood(A, predecessors, successors)
    println("Calcul des bornes F_bar et d_bar...")
    F_bar = Dict(j => instance_data.max_flying_time for j in fl_nodes)
    d_bar = Dict(j => d[j] for j in fl_nodes)
    for i in setdiff(V, fl_nodes)
        F_bar[i] = instance_data.max_flying_time
        d_bar[i] = 0
    end

    # Convergence pour F_bar (temps de vol maximum)
    converged = false
    while !converged
        converged = true
        # Set F_i ← min{F_i, max_{(i,j)∈A}{F_j - d_j}} for all i ∉ L_M^a
        for i in setdiff(fl_nodes, Set(L_a_M))  # i ∉ L_M^a
            successors_i = get(successors, i, [])  # Get flights j such that (i,j) ∈ A
            if !isempty(successors_i)
                old_Fi = F_bar[i]
                # Calculate max_{(i,j)∈A}{F_j - d_j}
                max_bound = maximum(F_bar[j] - d[j] for j in successors_i)
                F_bar[i] = min(F_bar[i], max_bound)
                if F_bar[i] != old_Fi
                    converged = false
                end
            end
        end
    end 
    # Convergence pour d_bar (durée minimum cumulée)
    converged = false
    iteration = 0
    
    while !converged
        converged = true
        iteration += 1
        #Set d_j ← max{d_j, min_{(i,j)∈A}{d_i + d_j}} for all j ∉ L_M
        for j in setdiff(fl_nodes, Set(L_M))  # j ∉ L_M
            predecessors_j = get(predecessors, j, [])  # Get flights i such that (i,j) ∈ A
            
            if !isempty(predecessors_j)
                old_dj = d_bar[j]
                # Calculate min_{(i,j)∈A}{d_i + d_j}
                min_bound = minimum(d_bar[i] + d[j] for i in predecessors_j)
                d_bar[j] = max(d_bar[j], min_bound)
                
                if d_bar[j] != old_dj
                    converged = false
                end
            end
        end
    end
    # Initialisation pour les nœuds non-vols
    #= for i in setdiff(V, fl_nodes)
        F_bar[i] = instance["maximum_flying_time"]
        d_bar[i] = 0
    end
        =#
    graph.arcs = A
    graph.nodes = V

    node_sets.fl_nodes = fl_nodes
    node_sets.ac_nodes = a_nodes
    node_sets.mtn_nodes = L_M  
    node_sets.mtn_nodes_stn = L_MS  

    arc_sets.arcs_K = A_K
    arc_sets.arcs_F = A_F
    arc_sets.arcs_S = A_S
    arc_sets.arcs_T = A_T
    arc_sets.arcs_M = A_M
    arc_sets.arcs_M_bar = A_M_bar

    fl_data.d = d
    fl_data.b = b
    fl_data.tk = tk
    fl_data.DT = DT
    fl_data.AT = AT

    other_data["M_FL_O"] = M_FL_O
    other_data["M_FL_D"] = M_FL_D   
    other_data["L_a_M"] = L_a_M
    other_data["a_day"] = a_day
    other_data["F_bar"] = F_bar
    other_data["d_bar"] = d_bar
    other_data["arc_reduc"] = arc_reduc
    other_data["node_reduc"] = node_reduc
    
    #= println(fl_nodes)
    println()
    println(instance_data.graph.node_sets.fl_nodes) =#

    return instance_data
end