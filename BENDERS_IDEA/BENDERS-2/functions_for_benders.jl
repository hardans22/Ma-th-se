function build_opti_cut(result, x, y, Z, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
    coupe = @build_constraint( sum(result.val_t1[(i,j)] * (max_flt * sum(x[(i,j),k] for k in a_nodes) - (max_flt - d[i]) * (1 - sum(y[j,k] for k in a_nodes))) for (i,j) in A_M) + 
            sum(result.val_t2[(i,j)] * (d[j] + (max_flt - d[i] - d[j]) * (1 - sum(x[(i,j),k] for k in a_nodes))) for (i,j) in A_theta2) + 
            sum(result.val_s[j] * (max_flt - (max_flt - d[j]) * sum(y[j,k] for k in a_nodes)) for j in L_M) + 
            sum(result.val_t3[(i,j)] * (d[j] - max_flt * (1 - sum(x[(i,j),k] for k in a_nodes))) for (i,j) in A_theta3) + 
            sum(result.val_t4[(i,j)] * (d[j] - max_flt * (1 - sum(x[(i,j),k] for k in a_nodes)) - max_flt * sum(y[j,k] for k in a_nodes)) for (i,j) in A_M) + 
            sum(result.val_p[k] * f[k] for k in a_nodes) + 
            sum(result.val_a1[j] * d[j] for j in V) + 
            sum(result.val_a2[j] * max_flt for j in V) <= Z)
    return coupe              
end

function build_fais_cut(result, x, y, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
    coupe = @build_constraint(sum(result.ray_t1[(i,j)] * (max_flt * sum(x[(i,j),k] for k in a_nodes) - (max_flt - d[i]) * (1 - sum(y[j,k] for k in a_nodes))) for (i,j) in A_M) + 
        sum(result.ray_t2[(i,j)] * (d[j] + (max_flt - d[i] - d[j]) * (1 - sum(x[(i,j),k] for k in a_nodes))) for (i,j) in A_theta2) + 
        sum(result.ray_s[j] * (max_flt - (max_flt - d[j]) * sum(y[j,k] for k in a_nodes)) for j in L_M) + 
        sum(result.ray_t3[(i,j)] * (d[j] - max_flt * (1 - sum(x[(i,j),k] for k in a_nodes))) for (i,j) in A_theta3) + 
        sum(result.ray_t4[(i,j)] * (d[j] - max_flt * (1 - sum(x[(i,j),k] for k in a_nodes)) - max_flt * sum(y[j,k] for k in a_nodes)) for (i,j) in A_M) + 
        sum(result.ray_p[k] * f[k] for k in a_nodes) + 
        sum(result.ray_a1[j] * d[j] for j in V) + 
        sum(result.ray_a2[j] * max_flt for j in V) <= 0)
    return coupe
end 

function compare_val(result, x_core, y_core, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
    val = sum(result.val_t1[(i,j)] * (max_flt * sum(x_core[(i,j,k)] for k in a_nodes) - (max_flt - d[i]) * (1 - sum(y_core[(j,k)] for k in a_nodes))) for (i,j) in A_M) + 
                    sum(result.val_t2[(i,j)] * (d[j] + (max_flt - d[i] - d[j]) * (1 - sum(x_core[(i,j,k)] for k in a_nodes))) for (i,j) in A_theta2) + 
                    sum(result.val_s[j] * (max_flt - (max_flt - d[j]) * sum(y_core[(j,k)] for k in a_nodes)) for j in L_M) + 
                    sum(result.val_t3[(i,j)] * (d[j] - max_flt * (1 - sum(x_core[(i,j,k)] for k in a_nodes))) for (i,j) in A_theta3) + 
                    sum(result.val_t4[(i,j)] * (d[j] - max_flt * (1 - sum(x_core[(i,j,k)] for k in a_nodes)) - max_flt * sum(y_core[(j,k)] for k in a_nodes)) for (i,j) in A_M) + 
                    sum(result.val_p[k] * f[k] for k in a_nodes) + 
                    sum(result.val_a1[j] * d[j] for j in V) + 
                    sum(result.val_a2[j] * max_flt for j in V)
    return val
end

function y_neighborhood(aircraft_paths, aircraft_mtn, y_val, acritik_set, A, L_M, A_M)
    y_list = []
    for ac in acritik_set
        ac_path = aircraft_paths[ac]
        current_mtn = aircraft_mtn[ac]
        if current_mtn != nothing
            current_idx = findfirst(==(current_mtn), ac_path)
            # Prendre seulement les nœuds de maintenance AVANT current_mtn dans le chemin
            #available_nodes = [ac_path[l] for l in 3:current_idx-1 if (ac_path[l-1], ac_path[l]) in A_M]
            available_nodes = [ac_path[l] for l in 3:length(ac_path) if (ac_path[l-1], ac_path[l]) in A_M && ac_path[l] != current_mtn]
            #= println("D'AUTRES OPPORTUNITÉS DE MAINTENANCE pour $ac")
            println(available_nodes)
             =#
            # Créer UNE solution pour CHAQUE nœud disponible
            for new_mtn in available_nodes
                y_new = copy(y_val)  # Copie de la solution actuelle
                y_new[current_mtn,ac] = 0  # Retirer l'ancienne maintenance
                y_new[new_mtn,ac] = 1      # Ajouter la nouvelle maintenance
                push!(y_list, y_new)  # Ajouter à la liste
            end
        end
    end
    return y_list
end

function maintenance_option(current_mtn, instance)
    L_M = instance["L_M"]
    DT = instance["DT"]
    end_horizon = instance["end_horizon_time"] 
    #return [j for j in L_M if DT[j] >= DT[current_mtn] && end_horizon - DT[current_mtn] < 1440 ]
    return [j for j in L_M if DT[j] >= DT[current_mtn]]
end 

function find_cover_set(aircraft_paths, y_val, acritik_set, succ_A, d, max_flt, L_M)
    all_cover_set = Vector{Tuple{Vector{String}, Int, String}}()
    for ac in acritik_set
        ac_path = aircraft_paths[ac]
        #= if aircraft_mtn[ac] != nothing
            continue    
        end =#
        #=  println("Aircraft critique pour cover set : $ac")
        println(ac_path)
         =#
        somme = 0
        ind = nothing 
        for indx in 1:length(ac_path)
            somme += d[ac_path[indx]]
            if somme > max_flt
                #= println(somme)
                println("Cover set pour l'avion $ac : ", (ac_path[1:indx], indx, ac))
                 =#
                ind = indx
                break
            end
        end
        if somme > max_flt
            #println(ind)
            cover_set = copy(ac_path[1:(ind-1)])
            key_node = ac_path[ind-1]
            somme -= d[ac_path[ind]]
             
            #println(somme)
            successors = get(succ_A, key_node, nothing)
            for j in successors
                if somme + d[j] > max_flt
                    push!(cover_set, j)
                    # Copier seulement si nécessaire
                    if sum(y_val[j,ac] for j in intersect(L_M, cover_set); init=0.0) >= 0.9
                        continue
                    end
                    push!(all_cover_set, (copy(cover_set), ind, ac))
                    pop!(cover_set)
                end
            end
        end 
    end 
    return all_cover_set
end


function print_information(spd_result, instance, A_theta2, A_theta3, eq_nodes)
    A_M, L_M = instance["A_M"], instance["L_M"]
    d, max_flt = instance["d"], instance["maximum_flying_time"] 
    a_nodes, f, V = instance["a_nodes"], instance["initial_flying_time"], instance["V"]
    
    if spd_result.status == MOI.DUAL_INFEASIBLE
        println("\n VARAIBLE THETA1")
        for a in A_M
            i,j = a
            if spd_result.ray_t1[a] != 0    
                println("i = $i, j = $j : $(spd_result.ray_t1[a])")
            end
        end

        println("\n VARAIBLE THETA2")
        for a in A_theta2
            i,j = a
            if spd_result.ray_t2[a] != 0    
                println("i = $i, j = $j : $(spd_result.ray_t2[a])")
            end
        end 

        println("\n VARAIBLE SIGMA")
        for j in L_M
            if spd_result.ray_s[j] != 0    
                println("j = $j : $(spd_result.ray_s[j])")
            end 
        end

        println("\n VARAIBLE THETA3")
        for a in A_theta3
            i,j = a
            if spd_result.ray_t3[a] != 0    
                println("i = $i, j = $j : $(spd_result.ray_t3[a])")
            end
        end

        println("\n VARAIBLE THETA4")        
        for a in A_M
            i,j = a
            if spd_result.ray_t4[a] != 0    
                println("i = $i, j = $j : $(spd_result.ray_t4[a])")
            end
        end

        println("\n VARAIBLE PI")
        for k in eq_nodes
            if spd_result.ray_p[k] != 0    
                println("k = $k : $(spd_result.ray_p[k])")
            end
        end

        println("\n VARAIBLE ALPHA1")
        for j in Set(V)
            if spd_result.ray_a1[j] != 0    
                println("j = $j : $(spd_result.ray_a1[j])")
            end
        end

        println("\n VARAIBLE ALPHA2")
        for j in Set(V)
            if spd_result.ray_a2[j] != 0    
                println("j = $j : $(spd_result.ray_a2[j])")
            end
        end
    else
        println("\n VARAIBLE THETA1")
        for a in A_M
            i,j = a
            if spd_result.val_t1[a] != 0    
                println("i = $i, j = $j : $(spd_result.val_t1[a])")
            end
        end

        println("\n VARAIBLE THETA2")
        for a in A_theta2
            i,j = a
            if spd_result.val_t2[a] != 0    
                println("i = $i, j = $j : $(spd_result.val_t2[a])")
            end
        end 

        println("\n VARAIBLE SIGMA")
        for j in L_M
            if spd_result.val_s[j] != 0    
                println("j = $j : $(spd_result.val_s[j])")
            end 
        end

        println("\n VARAIBLE THETA3")
        for a in A_theta3
            i,j = a
            if spd_result.val_t3[a] != 0    
                println("i = $i, j = $j : $(spd_result.val_t3[a])")
            end
        end

        println("\n VARAIBLE THETA4")        
        for a in A_M
            i,j = a
            if spd_result.val_t4[a] != 0    
                println("i = $i, j = $j : $(spd_result.val_t4[a])")
            end
        end

        println("\n VARAIBLE PI")
        for k in eq_nodes
            if spd_result.val_p[k] != 0    
                println("k = $k : $(spd_result.val_p[k])")
            end
        end

        println("\n VARAIBLE ALPHA1")
        for j in Set(V)
            if spd_result.val_a1[j] != 0    
                println("j = $j : $(spd_result.val_a1[j])")
            end
        end

        println("\n VARAIBLE ALPHA2")
        for j in Set(V)
            if spd_result.val_a2[j] != 0    
                println("j = $j : $(spd_result.val_a2[j])")
            end
        end
    end 
end

function build_aircraft_path(start_aircraft::String, succ::Dict{String, String}, d)
    """Construit le chemin complet d'un avion à partir du graphe des successeurs"""
    chemin = ["s", start_aircraft]
    current = start_aircraft
    cumul = d[start_aircraft]
    while haskey(succ, current) && succ[current] != "t"
        current = succ[current]
        cumul += d[current]
        push!(chemin, current)
    end
    
    # Ajouter le nœud terminal si accessible
    haskey(succ, current) && push!(chemin, succ[current])
    
    return chemin, cumul
end 

function build_all_paths(x_val, y_val, d, A, L_M, a_nodes)
    aircraft_paths = Dict{String, Vector{String}}()
    aircraft_mtn = Dict{String, Union{String, Nothing}}()
    aircraft_cumul = Dict{String, Float64}()
    s_rho = Dict(j => 0 for j in L_M)
    su = Dict()
    succ = Dict{String, String}()
    for (i, j) in Set(A)
        sum(x_val[(i,j,k)] for k in a_nodes) == 1.0 && (succ[i] = j)
    end
        
    for aircraft in a_nodes
        for j in L_M
            if y_val[j,aircraft] == 1.0
                aircraft_mtn[aircraft] = j
                break
            end
            aircraft_mtn[aircraft] = nothing

        end
        chemin, cumul = build_aircraft_path(aircraft, succ, d)
        aircraft_paths[aircraft] = chemin
        aircraft_cumul[aircraft] = cumul
    end
    return aircraft_paths, aircraft_mtn, aircraft_cumul
end

function find_aircraft_for_maintenance(maintenance_node::String, aircraft_paths::Dict{String, Vector{String}})
    """Trouve l'avion qui effectue la maintenance au nœud donné"""
    for (aircraft, path) in aircraft_paths
        maintenance_node in path && return aircraft
    end
    return "UNKNOWN"
end

function print_x_y(x_val, y_val, A, L_M, d, a_nodes)

    println("\nAFFICHAGE DE LA SOLUTION ")
    aircraft_paths = Dict{String, Vector{String}}()
    aircraft_mtn = Dict{String, String}()
    # Construire le graphe des successeurs en une seule passe
    for aircraft in a_nodes
        succ = Dict{String, String}()
        
        for (i, j) in A
            if x_val[(i,j,aircraft)] >= 0.9
                succ[i] = j
                #println("  Arc selected: $(i) -> $(j)")
            end
        end
    
        #println(succ)
        chemin, cumul = build_aircraft_path(aircraft, succ, d)
        println()
        #println(chemin)
        aircraft_paths[aircraft] = chemin
        cumul = 0.0
        path_with_cumul = String[]
        for i in 1:length(chemin)
            cumul += d[chemin[i]]
            push!(path_with_cumul, "$(chemin[i])_(d=$(d[chemin[i]]), Σ=$cumul)")
        end

        path_str = join(path_with_cumul, " ➡️  ")
        println("\n🛩️ Aircraft $aircraft path ($(length(chemin)-3) nodes): $path_str")
        println("   📊 Cumul total des d: $cumul")
        for j in L_M
            if y_val[j,aircraft] >= 0.9
                aircraft_mtn[aircraft] = j
            end
        end
        if haskey(aircraft_mtn, aircraft)
            println("\n🛠️ Maintenance at nodes: $(aircraft_mtn[aircraft]) by aircraft $aircraft")
        end
    
    end    
end


function find_cover_set(A, succ_A, x_val, y_val, a_nodes, d, max_flt, L_M)
    
    # 2. Pré-calculer L_M comme Set pour intersections rapides
    L_M_set = Set(L_M)
    
    # 3. Pré-allouer le vecteur de résultats avec estimation de taille
    all_cover_set = Vector{Tuple{Vector{String}, Int, String}}()
    sizehint!(all_cover_set, length(a_nodes) * 2)
    
    # 4. Réutiliser le buffer cover_set au lieu de créer/détruire
    cover_set = Vector{String}()
    sizehint!(cover_set, 20)  # Estimation taille typique
    
    @inbounds for aircraft in a_nodes
        succ = Dict{String, String}()
        sizehint!(succ, length(A))
        @inbounds for (i, j) in A
            if x_val[(i,j, aircraft)] >= 0.9
                succ[i] = j
            end
        end
        # Réinitialiser au lieu de recréer
        empty!(cover_set)
        push!(cover_set, "s", aircraft)
        
        current = aircraft
        temp = d[aircraft]
        len_cvs = 2
        
        # 5. Suivre le chemin avec get optimisé
        while true
            next = get(succ, current, nothing)
            
            # Break conditions combinées
            if isnothing(next) || next == "t"
                break
            end
            
            temp += d[next]
            push!(cover_set, next)
            len_cvs += 1
            current = next
            
            # Sortir dès qu'on dépasse
            if temp > max_flt
                if sum(y_val[j,aircraft] for j in intersect(L_M_set, cover_set); init=0.0) >= 0.9
                    break
                #= else
                    push!(all_cover_set, (copy(cover_set), len_cvs))  =#
                end
                 # 6. Traitement de la violation
                key_node = cover_set[end-1]
                temp -= d[cover_set[end]]
                pop!(cover_set)
                
                # 7. Obtenir successeurs une seule fois
                successors = get(succ_A, key_node, nothing)
                for j in successors
                    if temp + d[j] > max_flt
                        push!(cover_set, j)
                        # Copier seulement si nécessaire
                        push!(all_cover_set, (copy(cover_set), len_cvs, aircraft))
                        
                        pop!(cover_set)
                    end
                end 
                break
            end
        end
    end
    
    return all_cover_set
end


function find_cover_set_min(a_nodes, d, max_flt, succ_A)
    all_cv_set = []
    
    for start_node in a_nodes
        min_length = Inf
        node_results = []
        memo = Dict()  # Mémorisation
        
        function explore_paths(cv_set, somme, len_cvs, visited)
            # Clé pour mémorisation
            state_key = (cv_set[end], somme, len_cvs)
            if haskey(memo, state_key)
                return
            end
            memo[state_key] = true
            
            if somme > max_flt
                if len_cvs < min_length
                    min_length = len_cvs
                    empty!(node_results)
                    push!(node_results, (copy(cv_set), len_cvs, somme))
                elseif len_cvs == min_length
                    push!(node_results, (copy(cv_set), len_cvs, somme))
                end
                return
            end
            
            if len_cvs >= min_length
                return
            end
            
            current_node = cv_set[end]
            if current_node == "t"
                return
            end
            
            if haskey(succ_A, current_node)
                for next_node in succ_A[current_node]
                    if !(next_node in visited)
                        push!(cv_set, next_node)
                        push!(visited, next_node)
                        
                        explore_paths(cv_set, somme + d[next_node], len_cvs + 1, visited)
                        
                        pop!(cv_set)
                        pop!(visited)
                    end
                end
            end
        end
        
        cv_set = [start_node]
        visited = Set([start_node])
        explore_paths(cv_set, d[start_node], 1, visited)
        
        append!(all_cv_set, node_results)
    end
    
    return all_cv_set
end

function find_cover_set_2(a_nodes, d, max_flt, succ_A, nmax = 50)
    all_cv_set = []
    
    for start_node in a_nodes
        node_results = []
        count = 0
        memo = Dict()  # Mémorisation
        
        function explore_paths(cv_set, somme, len_cvs, visited)
            # Si on a déjà 100 chemins, arrêter
            if count >= 100
                return true
            end
            
            # Clé pour mémorisation (nœud actuel + somme)
            state_key = (cv_set[end], somme)
            if haskey(memo, state_key)
                return false
            end
            memo[state_key] = true
            
            # Si on dépasse max_flt, on enregistre
            if somme > max_flt
                push!(node_results, (copy(cv_set), len_cvs, somme))
                count += 1
                
                if count >= nmax
                    return true
                end
                return false  # Arrêter cette branche
            end
            
            # Obtenir les successeurs du dernier nœud
            current_node = cv_set[end]
            next_nodes = get(succ_A, current_node, nothing)
            
            if next_nodes !== nothing
                # Explorer tous les successeurs possibles
                for next_node in next_nodes
                    if next_node ∉ visited  # Éviter les cycles
                        push!(cv_set, next_node)
                        push!(visited, next_node)
                        
                        # Explorer et vérifier si on doit arrêter
                        should_stop = explore_paths(cv_set, somme + d[next_node], len_cvs + 1, visited)
                        
                        pop!(cv_set)
                        delete!(visited, next_node)
                        
                        # Si on a atteint 100 chemins, arrêter l'exploration
                        if should_stop
                            return true
                        end
                    end
                end
            end
            return false
        end
        
        cv_set = [start_node]
        visited = Set([start_node])
        explore_paths(cv_set, d[start_node], 1, visited)
        
        append!(all_cv_set, node_results)
    end
    
    return all_cv_set
end

function find_mtn_set(all_cover_set, A_M, d, max_flt)
    #all_mtn_set = Vector{NamedTuple{(:set_node, :len_set, :bnd_val), Tuple{Vector{String}, Int, Int}}}()
    best_mtn_set = (set_node = [], len_set = 0, bnd_val = 0)
    best_bnd_val = max_flt
    for (cover_set, len_cvs) in all_cover_set
        mtn_set = copy(cover_set)
        len_mtns = copy(len_cvs)
        end_node = mtn_set[end]
        while !((mtn_set[end-1], end_node) in A_M) && len_mtns > 2
            pop!(mtn_set)
            len_mtns -= 1
            end_node = mtn_set[end]
        end
        bnd_val = max_flt -  compute_fl_time(mtn_set, d) + d[end_node]
        #println("\nMAINTENANCE SET = $mtn_set")
        #println("BOUND = $bnd_val")
        #push!(all_mtn_set, (set_node = copy(mtn_set), len_set = len_mtns, bnd_val = bnd_val))
        if bnd_val < best_bnd_val
            best_bnd_val = bnd_val
            best_mtn_set = (set_node = copy(mtn_set), len_set = len_mtns, bnd_val = bnd_val)
        end
    end
    return best_mtn_set
end 

function compute_fl_time(set_nodes, d)
    total_time = 0
    for node in set_nodes
        total_time += d[node]
    end
    return total_time
end

