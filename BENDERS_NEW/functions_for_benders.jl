

function build_opti_cut(instance_data, sub_paths, x, Z, rho_val)
    fl_data = instance_data.fl_data
    max_flt = instance_data.max_flying_time
    f = fl_data.init_flying_time
    d = fl_data.d

    cut_list = []
    total_sum = AffExpr(0.0)  # Expression affine JuMP vide
    for (aircraft, (sub_path, len_sp)) in sub_paths
        #Première idée de coupe d'optimalité
        #= if len_sp >= 4
            sum_1 = sum((max_flt - d[sub_path[i]]) * (x[(sub_path[i], sub_path[i+1])] - 1) for i in 2:len_sp-2) - sum(d[sub_path[i+1]] * x[(sub_path[i], sub_path[i+1])] for i in 2:len_sp-2) - f[aircraft] + max_flt * x[(sub_path[end-1], sub_path[end])]  
        else 
            sum_1 = max_flt* x[(sub_path[2], sub_path[3])] - f[aircraft]
        end =#
        #Deuxième idée de coupe d'optimalité

        sum_1 = rho_val[sub_path[end]]*(sum(x[(sub_path[i],sub_path[i+1])] for i in 1:len_sp-1) - (len_sp-1) + 1)
        
        total_sum += sum_1
    end

    cut = @build_constraint(total_sum <= Z) 
    push!(cut_list, cut)
    return cut_list
end 

function build_feas_cut(infeas_subpaths, x)
    coupe_list = []
    for (_, (sub_path, len_sp)) in infeas_subpaths
        cut  = @build_constraint(sum(x[(sub_path[i],sub_path[i+1])] for  i in 1:len_sp-1) <= len_sp - 2)    
        push!(coupe_list, cut)
         
    end 
    return coupe_list
end 

function analyse_contrainte(c::ScalarConstraint)
    expr = c.func

    #Somme uniquement des coefficients des variables x
    somme_coeffs = sum(
        coeff for (var, coeff) in expr.terms 
        if startswith(name(var), "x"))

    rhs = c.set.upper

    println("Somme des coefficients de x : $somme_coeffs")
    println("Valeur droite (RHS)         : $rhs")
    println("Somme - RHS ?              : $(somme_coeffs - rhs)")

    return somme_coeffs, rhs
end

function infeasibility(model)

    # ⚠️ Toujours calculer le conflit avant
    compute_conflict!(model)

    iis_vars = VariableRef[]

    for (F, S) in list_of_constraint_types(model)
        for c in all_constraints(model, F, S)

            # Vérifie si la contrainte est dans le conflit
            status = MOI.get(model, MOI.ConstraintConflictStatus(), c)

            if status == MOI.IN_CONFLICT

                obj  = JuMP.constraint_object(c)
                func = obj.func

                # ==========================
                # Cas expression affine
                # ==========================
                if func isa AffExpr
                    for term in func.terms
                        var = term.first

                        # 🔹 garder uniquement les variables u
                        if startswith(name(var), "u")
                            push!(iis_vars, var)
                        end
                    end
                end

            end
        end
    end

    return unique(iis_vars)
end

function extract_u_indices(iis_vars)

    indices = String[]

    for var in iis_vars
        idx = JuMP.name(var)      # ex: "u[N121]"

        # enlever le "u[" au début et "]" à la fin
        clean = replace(idx, "u[" => "")
        clean = replace(clean, "]" => "")

        push!(indices, clean)
    end

    return indices
end

function check_constraints(constraints_list, solution; tol=1e-6)

    for con in constraints_list
        println()
        println(con)
        func = con.func
        set  = con.set

        # calcul de la valeur de l'expression
        lhs = func.constant

        for term in func.terms
            println(term.variable)
            lhs += term.coefficient * solution[term.variable]
        end

        if set isa MOI.LessThan
            if lhs > set.upper + tol
                println("❌ Contrainte violée :", con)
                return false
            end

        elseif set isa MOI.GreaterThan
            if lhs < set.lower - tol
                println("❌ Contrainte violée :", con)
                return false
            end

        elseif set isa MOI.EqualTo
            if abs(lhs - set.value) > tol
                println("❌ Contrainte violée :", con)
                return false
            end
        end
    end

    return true
end

function compare_sp_sol(result, result_1, instance_data)
    graph = instance_data.graph
    node_sets = graph.node_sets
    
    L_M = node_sets.mtn_nodes
    V = graph.nodes

    y_val, u_val, rho_val = result.y, result.u, result.rho
    y_val_1, u_val_1, rho_val_1 = result_1.y, result_1.u, result_1.rho

    for j in L_M
        if abs(y_val[j] - y_val_1[j]) >= 1e-6 && abs(rho_val[j] - rho_val_1[j]) >= 1e-6
            println("IL S'AGIT DE Y ET RHO")
            return false
        end
    end
    for i in V
        if abs(u_val[i] - u_val_1[i]) >= 1e-6
            println("IL S'AGIT DE U")
            println("i = ", i)
            println("u_val[i] = ", u_val[i])
            println("u_val_1[i] = ", u_val_1[i])
            return false
        end 
    end 
    return true 
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
                println("j = $j : $(spd_sub_path[end-1]result.ray_a2[j])")
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

function build_aircraft_path(start_aircraft::String, succ::Dict{String, String})
    """Construit le chemin complet d'un aircraft à partir du graphe des successeurs"""
    chemin = ["s", start_aircraft]
    current = start_aircraft
    
    while haskey(succ, current) && succ[current] != "t"
        current = succ[current]
        push!(chemin, current)
    end
    
    # Ajouter le nœud terminal si accessible
    haskey(succ, current) && push!(chemin, succ[current])
    
    return chemin
end 

function build_all_paths(x_val, A)
    """Construit les chemins de tous les avions à partir des variables x"""
    succ = Dict{String, String}()
    for (i, j) in Set(A)
        x_val[(i, j)] >= 0.9 && (succ[i] = j)
    end
    
    aircraft_paths = Dict{String, Vector{String}}()
    
    # Identifier les avions et construire leurs chemins
    aircraft_starts = [(i, j) for (i, j) in A if i == "s" && x_val[(i, j)] >= 0.9]
    
    for (_, avion) in aircraft_starts
        chemin = build_aircraft_path(avion, succ)
        aircraft_paths[avion] = chemin
    end
    
    return aircraft_paths
end

function find_aircraft_for_maintenance(maintenance_node::String, aircraft_paths::Dict{String, Vector{String}})
    """Trouve l'aircraft qui effectue la maintenance au nœud donné"""
    for (aircraft, path) in aircraft_paths
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
    aircraft_paths = Dict{String, Vector{String}}()
    longest_aircraft = ""
    longest_path = String[]
    longest_length = 0
    
    # Identifier les aircrafts et construire leurs chemins
    #aircraft_starts = [("s", j) for j in a_nodes if x_val[("s", j)] >= 0.9]
    
    println("✈️ Paths for each aircraft:")
    aircraft_mtn = Dict{String, Vector{String}}()
    for aircraft in H_aircraft
        chemin = build_aircraft_path(aircraft, succ)
        aircraft_paths[aircraft] = chemin
        
        # Vérifier si c'est le chemin le plus long
        if length(chemin) > longest_length
            longest_length = calculate_path_distance(chemin, d)
            longest_path = chemin
            longest_aircraft = aircraft
        end
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


