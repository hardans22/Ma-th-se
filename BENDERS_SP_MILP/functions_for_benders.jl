include("../structures.jl")

function build_opti_cut(sub_paths, x, Z, rho_val)
    cut_list = []
    sum_1 = AffExpr(0.0) 
    for (aircraft, (sub_path, len_sp, mtn_node)) in sub_paths
        if rho_val[mtn_node] > 0
            sum_1 = rho_val[mtn_node]*(1 - sum((1 - x[(sub_path[i],sub_path[i+1])]) for i in 1:len_sp-1))
            cut = @build_constraint(sum_1 <= Z[aircraft]) 
            push!(cut_list, cut)
        end
    end
    return cut_list
end 

function build_opti_cut_agg(sub_paths, x, Z, obj_sp)
    summ = AffExpr(0.0) 
    big_set = []
    big_set_len = 0.0 
    for (aircraft, (sub_path, len_sp)) in sub_paths
        if len_sp != 0
            summ += 1 - sum((1 - x[(sub_path[i],sub_path[i+1])]) for i in 1:len_sp-1)
        end
    end
    cut = @build_constraint(obj_sp*(summ) <= Z) 
    return cut  
end 


function build_feas_cut(infeas_subpaths, x)
    cut_list = []
    for (_, (sub_path, len_sp)) in infeas_subpaths
        cut  = @build_constraint(sum(x[(sub_path[i],sub_path[i+1])] for  i in 1:len_sp-1) <= len_sp - 2)    
        push!(cut_list, cut)
    end 
    return cut_list
end 


function irreductible_set(sub_paths, instance_data, obj_sp, d)
    aircraft_paths = Dict(aircraft => sub_path for (aircraft, (sub_path, len_sp, mtn_node)) in sub_paths)
    a_nodes = keys(aircraft_paths)
    #= result = solve_sp_pd(aircraft_paths, instance_data, a_nodes)
    println()
    println("Obj du SP = $obj_sp")
    println("Obj du SP réduit = $(result.obj)")
    println(print_u(result.u, aircraft_paths))
    =#
    result = solve_sp_pd(aircraft_paths, instance_data, a_nodes)
    println(print_u(result.u, aircraft_paths))



    for (key, sub_path) in aircraft_paths
        new_sub_path = sub_path[1:end-1]
        aircraft_paths[key] = new_sub_path
        result = solve_sp_pd(aircraft_paths, instance_data, a_nodes)
        #= println("Obj du SP = $obj_sp")
        println("Obj après un élément supprimé = $(result.obj)")
         =#
        if result.obj != obj_sp
            #= println("YYYYYEEEEEEEEEEEEEEESSSSSSSSSS")
            println(print_u(result.u, aircraft_paths)) =#
            break
        else
            #= for key in keys(sub_paths)
                println()
                println(key)
                sub_path = sub_paths[key][1]
                println("Somme du sous chemin = ", sum(d[i] for i in sub_path))
            end
             =#
            println("APRÈS SUPPRESSION")
            #= for key in keys(aircraft_paths)
                println()
                println(key)
                sub_path = aircraft_paths[key]
                println("Somme du sous chemin = ", sum(d[i] for i in sub_path))
            end
             =#

            println("NNNNNNNNOOOOOOOOOOOOOOOOOOOOONNNNNNNNNNN")
            println("Obj du SP = $obj_sp")
            println("Obj après un élément supprimé = $(result.obj)")
            println(print_u(result.u, aircraft_paths))
            continue
        end
    end 
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


