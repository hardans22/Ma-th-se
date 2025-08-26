using JSON

# ====================================================================
# FONCTIONS UTILITAIRES POUR LA MANIPULATION DE GRAPHES
# ====================================================================

"""Retourne tous les arcs sortants du nœud i"""
delta_plus(i, A) = [(u, v) for (u, v) in A if u == i]

"""Retourne tous les arcs entrants vers le nœud i"""
delta_minus(i, A) = [(u, v) for (u, v) in A if v == i]

"""Retourne la liste des prédécesseurs du nœud i"""
predecessor(i, A) = [x[1] for x in A if x[2] == i]

"""Retourne la liste des successeurs du nœud i"""
successor(i, A) = [x[2] for x in A if x[1] == i]

"""
Met à jour les ensembles d'arcs en remplaçant 'arc' par 'new_arc'
Utilise des Sets temporaires pour des recherches O(1)
"""
function arc_set_update(arc, new_arc, A, A_S, A_T, A_M)
    push!(A, new_arc)
    
    # Conversion en Sets pour recherches rapides
    sets = (Set(A_S), Set(A_T))
    arrays = (A_S, A_T)
    
    # Mise à jour dans le bon ensemble
    for (s, arr) in zip(sets, arrays)
        if arc in s
            #println(arc)
            push!(arr, new_arc)
            filter!(x -> x != arc, arr)
            break
        end
    end
    
    # Traitement spécial pour A_M
    if arc in Set(A_M)
        push!(A_M, new_arc)
        filter!(x -> x != arc, A_M)
    end
    
    return A, A_S, A_T, A_M
end

"""
Supprime un arc des différents ensembles d'arcs
"""
function rm_arcs(arc, A_K, A_S, A_F, A_T, A_M)
    sets = (Set(A_K), Set(A_S), Set(A_F), Set(A_T))
    arrays = (A_K, A_S, A_F, A_T)
    
    # Suppression dans le bon ensemble
    for (s, arr) in zip(sets, arrays)
        if arc in s
            filter!(x -> x != arc, arr)
            break
        end
    end
    
    # Traitement spécial pour A_M
    arc in Set(A_M) && filter!(x -> x != arc, A_M)
    
    return A_K, A_S, A_F, A_T, A_M
end

"""
Met à jour un ensemble de nœuds lors de la fusion de i et j en ij
"""
function nodes_set_update(i, j, ij, nodes_set)
    i in nodes_set && (nodes_set = setdiff(nodes_set, [i]); push!(nodes_set, ij))
    j in nodes_set && (nodes_set = setdiff(nodes_set, [j]))
    return nodes_set
end

"""
Met à jour les données (durées, temps) lors de la fusion de deux nœuds
"""
function data_update(i, j, ij, d, DT, AT, nodes_rm)
    # Fusion des données
    d[ij] = d[i] + d[j]
    DT[ij] = DT[i]
    AT[ij] = AT[j]
    # Suppression des anciens nœuds (sauf source et puits)
    for (node, name) in [(i, "s"), (j, "t")]
        if node != name
            push!(nodes_rm, node)
            for dict in [d, DT, AT]
                pop!(dict, node, nothing)
            end
        end
    end
    return d, DT, AT, nodes_rm
end


# ====================================================================
# FONCTION PRINCIPALE DE CONSTRUCTION DU GRAPHE
# ====================================================================

"""
Construit le graphe de planification des vols avec preprocessing optionnel
- file: fichier JSON contenant l'instance
- preprocess: si true, applique la réduction de graphe
"""
function build_graph(file, preprocess)
    instance = JSON.parsefile(file)
    
    # ================================================================
    # EXTRACTION DES DONNÉES D'ENTRÉE
    # ================================================================
    
    nbr_FL, nbr_K = instance["number_of_flight_legs"], instance["number_of_aircrafts"]
    f_legs, aircrafts = instance["flight_legs"], instance["aircrafts"]
    mtn_stations = Set(instance["maintenance_stations"])  # Set pour recherche O(1)
    init_airport = instance["initial_airport_aircraft"]
    TRT, mtn_time = instance["turn_around_time"], instance["maintenance_time"]
    initial_flying_time = instance["initial_flying_time"]
    initial_takeoff, initial_flying_day = instance["initial_takeoff"], instance["initial_flying_day"]
    end_h_time, nbr_TP = instance["end_horizon_time"], instance["nbr_TP"]
    
    # ================================================================
    # CONSTRUCTION DES NŒUDS
    # ================================================================
    
    # Nœuds de vols : "origine_destination_départ_arrivée"
    fl_nodes = [f_legs[i][1]*"_"*f_legs[i][2]*"_"*string(f_legs[i][3])*"_"*string(f_legs[i][4]) 
                for i in 1:nbr_FL]
    
    # Tous les nœuds
    a_nodes = aircrafts
    st_nodes = ["s", "t"]
    V = vcat(fl_nodes, a_nodes, st_nodes)
    V_wt_st = vcat(a_nodes, fl_nodes)
    
    # ================================================================
    # INITIALISATION DES STRUCTURES DE DONNÉES
    # ================================================================
    
    # Dictionnaires de données des vols
    flight_legs = Dict(fl_nodes[i] => [f_legs[i][1], f_legs[i][2]] for i in 1:nbr_FL)
    DT = Dict(fl_nodes[i] => f_legs[i][3] for i in 1:nbr_FL)  # Departure Time
    AT = Dict(fl_nodes[i] => f_legs[i][4] for i in 1:nbr_FL)  # Arrival Time
    d = Dict(fl_nodes[i] => f_legs[i][5] for i in 1:nbr_FL)   # Duration
    
    # Initialisation des avions et nœuds source/puits
    A_S = [(st_nodes[1], k) for k in a_nodes]  # Arcs source -> avions
    
    for k in a_nodes
        d[k], DT[k], AT[k] = initial_flying_time[k], 0, 0
    end
    
    for i in st_nodes
        d[i], DT[i], AT[i] = 0, 0, 0
    end
    
    # ================================================================
    # CONSTRUCTION DES ARCS
    # ================================================================
    
    A_K, A_F, A_T = [], [], []  # Arcs avions->vols, vols->vols, vols->puits
    a = Dict()  # Matrice d'assignation temporelle
    
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
        for j in fl_nodes
            if i != j && dest == flight_legs[j][1]
                time_diff = DT[j] - at_i
                if TRT <= time_diff <= 1440
                    push!(A_F, (i, j))
                end
            end
        end
        
        # A_T: connexions vol -> puits
        (end_h_time - at_i <= 1440) && push!(A_T, (i, st_nodes[2]))
        
        # Matrice d'assignation par période temporelle
        day = ceil(Int, dt_i / 1440)
        for t in 1:nbr_TP
            a[(i, t)] = (day == t) ? 1 : 0
        end
    end
    
    # ================================================================
    # ARCS DE MAINTENANCE ET VARIABLES BINAIRES
    # ================================================================
    
    A = vcat(A_S, A_K, A_F, A_T)
    old_arc_nbr = length(A)
    A_M, L_M = [], []  # Arcs et nœuds de maintenance
    b = Dict()  # Variable binaire jour
    
    println("Identification des opportunités de maintenance...")
    for x in A
        i, j = x
        # Détermination de la station d'arrivée
        arrival_station = if i in fl_nodes
            flight_legs[i][2]
        elseif i in a_nodes
            init_airport[i]
        else
            ""
        end
        
        # Arc de maintenance possible ?
        if j in fl_nodes && arrival_station in mtn_stations && DT[j] - AT[i] >= mtn_time
            push!(A_M, x)
            push!(L_M, j)
        end
        
        # Variable binaire pour changement de jour
        b[x] = (floor(Int, DT[j]/1440) - floor(Int, DT[i]/1440) == 1) ? 1 : 0
    end
    
    L_M = unique(L_M)
    A_M_bar = setdiff(A, A_M)
    
    # Regroupement des vols par station de maintenance
    L_MS = Dict(ms => [i for i in L_M if flight_legs[i][1] == ms] for ms in mtn_stations)
    
    println("A AVANT : $(length(A)) arcs")
    
    # ================================================================
    # PREPROCESSING : RÉDUCTION DU GRAPHE
    # ================================================================
    
    if preprocess
        println("Début du preprocessing...")
        queue = copy(A_M_bar)
        visited = Set()
        
        while !isempty(queue)
            arc = popfirst!(queue)
            arc in visited && continue
            push!(visited, arc)
            
            i, j = arc
            arcs_rm, nodes_rm, nodes_new = [], [], []
            ij = i * "_" * j
            
            # --------------------------------------------------------
            # CAS 1: i est l'unique prédécesseur de j
            # --------------------------------------------------------
            if !(arc in Set(A_S)) && predecessor(j, A) == [i]
                outarc_i, inarc_i, outarc_j = delta_plus(i, A), delta_minus(i, A), delta_plus(j, A)
                push!(nodes_new, ij)
                
                # Mise à jour des ensembles de nœuds
                if i in fl_nodes
                    fl_nodes = setdiff(fl_nodes, [i])
                    push!(fl_nodes, ij)
                elseif i in a_nodes
                    a_nodes = setdiff(a_nodes, [i])
                    push!(a_nodes, ij)
                    
                    # Transfert des propriétés d'avion
                    initial_flying_time[ij] = initial_flying_time[i] + d[j]
                    initial_takeoff[ij] = initial_takeoff[i] + 1
                    initial_flying_day[ij] = initial_flying_day[i]
                    
                    for dict in [initial_flying_time, initial_takeoff, initial_flying_day]
                        filter!(kv -> kv[1] != i, dict)
                    end
                end
                
                fl_nodes = setdiff(fl_nodes, [j])
                
                # Redirection des arcs
                for inarc in inarc_i
                    new_arc = (inarc[1], ij)
                    A, A_S, A_T, A_M = arc_set_update(inarc, new_arc, A, A_S, A_T, A_M)
                    push!(queue, new_arc)
                end
                
                for outarc in outarc_j
                    new_arc = (ij, outarc[2])
                    A, A_S, A_T, A_M = arc_set_update(outarc, new_arc, A, A_S, A_T, A_M)
                    push!(queue, new_arc)
                end
                
                # Suppression des anciens arcs
                for outarc in outarc_i
                    rm_arcs(outarc, A_K, A_S, A_F, A_T, A_M)
                end
                
                append!(arcs_rm, vcat(outarc_i, inarc_i, outarc_j))
                d, DT, AT, nodes_rm = data_update(i, j, ij, d, DT, AT, nodes_rm)
                
                # Mise à jour des structures de maintenance
                L_M = nodes_set_update(i, j, ij, L_M)
                for ms in mtn_stations
                    L_MS[ms] = nodes_set_update(i, j, ij, L_MS[ms])
                end
            end
            
            # Nettoyage après cas 1
            A = setdiff(A, arcs_rm)
            A_M_bar = setdiff(A, A_M)
            V = setdiff(V, nodes_rm)
            append!(V, nodes_new)
            unique!(V)
            
            # --------------------------------------------------------
            # CAS 2: j est l'unique successeur de i
            # --------------------------------------------------------
            arcs_rm, nodes_rm, nodes_new = [], [], []
            
            if !(arc in Set(A_T)) && successor(i, A) == [j]
                inarc_j, inarc_i, outarc_j = delta_minus(j, A), delta_minus(i, A), delta_plus(j, A)
                push!(nodes_new, ij)
                
                # Logique similaire au cas 1 (mise à jour des nœuds et arcs)
                if i in fl_nodes
                    fl_nodes = setdiff(fl_nodes, [i])
                    push!(fl_nodes, ij)
                elseif i in a_nodes
                    a_nodes = setdiff(a_nodes, [i])
                    push!(a_nodes, ij)
                    
                    initial_flying_time[ij] = initial_flying_time[i] + d[j]
                    initial_takeoff[ij] = initial_takeoff[i] + 1
                    initial_flying_day[ij] = initial_flying_day[i]
                    
                    for dict in [initial_flying_time, initial_takeoff, initial_flying_day]
                        filter!(kv -> kv[1] != i, dict)
                    end
                end
                
                fl_nodes = setdiff(fl_nodes, [j])
                
                # Redirection des arcs
                for inarc in inarc_i
                    new_arc = (inarc[1], ij)
                    A, A_S, A_T, A_M = arc_set_update(inarc, new_arc, A, A_S, A_T, A_M)
                    push!(queue, new_arc)
                end
                
                for outarc in outarc_j
                    new_arc = (ij, outarc[2])
                    A, A_S, A_T, A_M = arc_set_update(outarc, new_arc, A, A_S, A_T, A_M)
                    push!(queue, new_arc)
                end
                
                for inarc in inarc_j
                    rm_arcs(inarc, A_K, A_S, A_F, A_T, A_M)
                end
                
                append!(arcs_rm, vcat(inarc_j, inarc_i, outarc_j))
                d, DT, AT, nodes_rm = data_update(i, j, ij, d, DT, AT, nodes_rm)
                
                L_M = nodes_set_update(i, j, ij, L_M)
                for ms in mtn_stations
                    L_MS[ms] = nodes_set_update(i, j, ij, L_MS[ms])
                end
            end
            
            # Nettoyage après cas 2
            A = setdiff(A, arcs_rm)
            A_M_bar = setdiff(A, A_M)
            V = setdiff(V, nodes_rm)
            append!(V, nodes_new)
            unique!(V)
        end
        
        # ============================================================
        # MISE À JOUR POST-PREPROCESSING
        # ============================================================
        
        V_wt_st = vcat(fl_nodes, a_nodes)
        unique!(V_wt_st)
        println("A APRÈS : $(length(A)) arcs")
        
       
        
        A_M_bar = setdiff(A, A_M)
        
        # Recalcul des variables binaires et structures
        b = Dict()
        L_a_M = []
        A_K = []
        A_F = []
        for arc in Set(A)
            i, j = arc[1], arc[2] 
            b[arc] = (floor(Int, DT[j]/1440) - floor(Int, DT[i]/1440) == 1) ? 1 : 0
            j in L_M && push!(L_a_M, i)
            if i in Set(a_nodes) 
                push!(A_K, arc)
            end
            if i in Set(fl_nodes) && j in Set(fl_nodes)
                push!(A_F, arc)
            end  
        end
        #println(A_K)
        unique!(L_a_M)
        
        # Recalcul de la matrice d'assignation temporelle
        a = Dict()
        for i in fl_nodes
            day = ceil(Int, DT[i] / 1440)
            for t in 1:nbr_TP
                a[(i, t)] = (day == t) ? 1 : 0
            end
        end
        
        # ============================================================
        # CALCULS DE CONVERGENCE : F_bar et d_bar
        # ============================================================
        
        println("Calcul des bornes F_bar et d_bar...")
        F_bar = Dict(j => instance["maximum_flying_time"] for j in fl_nodes)
        d_bar = Dict(j => d[j] for j in fl_nodes)
        L_a_M_set = Set(L_a_M)
        
        # Convergence pour F_bar (temps de vol maximum)
        converged = false
        while !converged
            converged = true
            # Set F_i ← min{F_i, max_{(i,j)∈A}{F_j - d_j}} for all i ∉ L_M^a
            for i in setdiff(fl_nodes, L_a_M)  # i ∉ L_M^a
                successors_i = successor(i, A_F)  # Get flights j such that (i,j) ∈ A
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
            for j in setdiff(fl_nodes, L_M)  # j ∉ L_M
                predecessors_j = predecessor(j, A_F)  # Get flights i such that (i,j) ∈ A
                
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
        for i in setdiff(V, fl_nodes)
            F_bar[i] = instance["maximum_flying_time"]
            d_bar[i] = 0
        end

        instance["F_bar"], instance["d_bar"] = F_bar, d_bar
    end
    F_bar = Dict(j => instance["maximum_flying_time"] for j in V)
    d_bar = Dict(j => d[j] for j in V)
    instance["F_bar"], instance["d_bar"] = F_bar, d_bar

    # ================================================================
    # FINALISATION ET RETOUR
    # ================================================================
    
    new_arc_nbr = length(A)
    arc_reduc = old_arc_nbr - new_arc_nbr
    
    # Mise à jour de l'instance avec toutes les structures calculées
    merge!(instance, Dict(
        "DT" => DT, "AT" => AT, "d" => d, "b" => b,
        "A_S" => A_S, "A_F" => A_F, "A_K" => A_K, "A_T" => A_T,
        "A" => A, "A_M" => A_M, "A_M_bar" => A_M_bar, "L_M" => L_M,
        "L_MS" => L_MS, "V" => V, "V_wt_st" => V_wt_st, "a" => a,
        "a_nodes" => a_nodes, "fl_nodes" => fl_nodes,
        "initial_flying_time" => initial_flying_time,
        "initial_takeoff" => initial_takeoff,
        "initial_flying_day" => initial_flying_day,
        "arc_reduc" => arc_reduc
    ))
    
    println("Construction terminée. Réduction: $arc_reduc arcs")
    return instance
end