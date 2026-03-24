include("functions_for_benders.jl")
include("spdual_functions.jl")

function benders_decomp(env, instance_file, output_file, nbr_thread, silent, time_limit)
    instance = build_graph(instance_file, false)
    nbr_K = instance["number_of_aircrafts"]
    max_flt = instance["maximum_flying_time"]
    d = instance["d"]
    A_S, A_T = instance["A_S"], instance["A_T"]
    A_F, A_K = instance["A_F"], instance["A_K"]
    A_M, A, L_M, A_M_bar = instance["A_M"], instance["A"], instance["L_M"], instance["A_M_bar"]
    V_wt_st, V, a_nodes = instance["V_wt_st"], instance["V"], instance["a_nodes"]
    f = instance["initial_flying_time"]
    fl_nodes = instance["fl_nodes"]
    mtn_stations, ms_capacity = instance["maintenance_stations"], instance["mtn_station_capacity"]
    MS, L_MS = instance["maintenance_stations"], instance["L_MS"]
    ac_critique  = instance["ac_critique"]
    acritik_set = instance["acritik_set"]
    
    
    pred_A, succ_A = Dict(), Dict()
    pred_A, succ_A = update_neighborhood(A, pred_A, succ_A)
    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(A_M, pred_A_M, succ_A_M)
    pred_A_M_bar, succ_A_M_bar = Dict(), Dict()
    pred_A_M_bar, succ_A_M_bar = update_neighborhood(A_M_bar, pred_A_M_bar, succ_A_M_bar)
    
    #= @time begin
        all_cv_set_min = find_cover_set_min(fl_nodes, acritik_set, d, max_flt, pred_A, succ_A)
    end 
    println("NOMBRE DE COVER SET = ", length(all_cv_set_min)) 

    for cv_set in all_cv_set_min
        println()
        println(cv_set)
    end  =#

    
    # ===================== Master problem =====================
    
    master_model = Model(() -> Gurobi.Optimizer(env))
    GRBsetintparam(env, "OutputFlag", silent ? 0 : 1)
    set_optimizer_attribute(master_model, "Threads", nbr_thread)
    set_optimizer_attribute(master_model, "TimeLimit", time_limit) 
    set_optimizer_attribute(master_model, "LazyConstraints", 1)
    set_optimizer_attribute(master_model, "Cuts", 2)
    set_optimizer_attribute(master_model, "IntFeasTol", 1e-9)
    set_optimizer_attribute(master_model, "LogFile", output_file)

    #set_optimizer_attribute(master_model, "MIPGap", 1e-9)
    #set_optimizer_attribute(master_model, "DisplayInterval", 1)           # Aggressive
    # Variables
    @variable(master_model, x[a in A, k in a_nodes], Bin)
    @variable(master_model, y[j in L_M, k in a_nodes], Bin)
    @variable(master_model, Z >= 0)
    master_model[:x] = x
    master_model[:y] = y
    #= for k in acritik_set, a in A
        MOI.set(master_model, Gurobi.VariableAttribute("BranchPriority"), x[a,k], 100)
    end
     =##= for a in A
        println(a)
    end =#
    A_M_bar_without_t = [a for a in A_M_bar if a[2] != "t"]
    NH_aircraft = setdiff(a_nodes, acritik_set)
    @constraint(master_model, mp_c1[k in a_nodes], x[("s",k),k] == 1)
    @constraint(master_model, mp_c2[k in a_nodes], sum(x[a,k] for a in A_T) == 1)
    @constraint(master_model, mp_c3[i in fl_nodes], sum(x[(i,j),k] for j in get(succ_A, i, []), k in a_nodes) == 1)
    @constraint(master_model, mp_c4[i in V_wt_st, k in a_nodes], sum(x[(j,i),k] for j in get(pred_A, i, [])) == sum(x[(i,j),k] for j in get(succ_A, i, [])))
    @constraint(master_model, mp_c5[j in L_M, k in a_nodes], y[j,k] <= sum(x[(i,j),k] for i in get(pred_A_M, j, [])))

    @constraint(master_model, mp_c6[j in L_M, k in NH_aircraft], sum(y[j,k] for j in L_M) == 0)
    #@constraint(master_model, sum(y[j,k] for j in L_M, k in a_nodes) <= ac_critique)   #Pas necessaire
    @constraint(master_model, mp_c7[k in acritik_set], sum(y[j,k] for j in L_M) <= 1)
    @constraint(master_model, mp_c8[k in NH_aircraft], x[("s",k),k] + sum(d[j]*x[(k,j),k] for j in get(succ_A, k, [])) + sum(d[j]*x[(i,j),k] for (i,j) in A_F) <= max_flt)
    @constraint(master_model, mp_c9[k in a_nodes], sum(x[(k,j),k] for j in get(succ_A, k, [])) == 1)
    #println(mp_c8)

    #= @constraint(master_model, mp_c10[j in L_M], sum(y[j,k] for k in a_nodes) <= 1)
    @constraint(master_model, sum(y[j,k] for j in L_M for k in acritik_set) >= 1)
     =#
    
    @constraint(master_model, mp_c11[k in acritik_set], d[k]*x[("s",k),k] + sum(d[j]*x[(k,j),k] for j in get(succ_A, k, [])) + sum(d[j]*x[(i,j),k] for (i,j) in A_F) <= max_flt*(1+sum(y[j,k] for j in L_M)))
    #@constraint(master_model, Z <= 5)
    #write_to_file(master_model, "master_model.lp")
    
    # Contrainte de brisure de symétrie seulement sur les avions non-critiques
    
    #Valid inequalities
    #@constraint(master_model, c40[i in V], sum(d[j]*x[(i,j)] for j in get(succ_A, i, [])) <= max_flt*(1 + sum(y[j] for j in L_M)))
    #@constraint(master_model,sum(sum(d[j]*x[(i,j)] for j in get(succ_A, i, [])) for i in V) <= max_flt*nbr_K*(1 + sum(y[j] for j in L_M)))
    
    #= @constraint(master_model, c32[j in V_wt_st], (length(get(pred_A_M_bar, j, []))-1) + sum(x[(i, j)] for i in get(pred_A_M_bar, j, [])) >= 1)
    @constraint(master_model, c33[i in V_wt_st], (length(get(succ_A_M_bar, i, []))-1) + sum(x[(i, j)] for j in get(succ_A_M_bar, i, [])) >= 1)
     =# 
    #@objective(master_model, Min, sum(y[j] for j in L_M) + Z)
    @objective(master_model, Min, Z)
    
    #= prep1_time = @elapsed begin
          
        all_cover_set = find_cover_set_2(acritik_set, d, max_flt, succ_A, 100)

        println("NOMBRE DE COVER SET = ", length(all_cover_set)) 
        added_cover_sets = []
        for element in all_cover_set 
            #println(element)
            cover_set = element[1]
            len_cvs = element[2]
            for ac in acritik_set
                if cover_set in added_cover_sets
                    #println("DÉJÀ TROUVÉE")
                    continue
                else 
                    arc_cover_set = [(cover_set[i],cover_set[i+1]) for i in 1:len_cvs-1]
                    mtn_possibilities = intersect(A_M, arc_cover_set)
                    j_poss = unique([arc[2] for arc in mtn_possibilities])
                    @constraint(master_model, sum(x[a,ac] for a in arc_cover_set) <= (len_cvs - 2)*(1 + sum(y[j,ac] for j in j_poss)))
                    push!(added_cover_sets, cover_set)
                end
            end
        end
    end

    prep2_time = @elapsed begin
        set_optimizer_attribute(master_model, "PoolSearchMode", 2)  # Mode exhaustif
        set_optimizer_attribute(master_model, "PoolSolutions",50)  # Nombre max de solutions à stocker
        set_optimizer_attribute(master_model, "PoolGap", 0.1)  # Solutions jusqu'à 20% de l'optimal

        # Pour forcer plus de diversité
        #set_optimizer_attribute(master_model, "SolutionLimit", 10)

        optimize!(master_model)
        n_solutions = result_count(master_model)
        println("\nNOMBRE DE SOLUTION OPTIMALES : ", n_solutions)
    
        # ÉTAPE 1 : D'abord, récupérer TOUTES les solutions sans modifier le modèle
        all_solutions = []
        for i in 1:n_solutions
            x_values = value.(x; result = i)
            y_values = value.(y; result = i)
            push!(all_solutions, (x = x_values, y = y_values))
        end

        # ÉTAPE 2 : Maintenant, traiter les solutions et modifier le modèle
        compt_cst = 0
        added_cover_sets = []
        for (i, sol) in enumerate(all_solutions)
            #println("\nSolution $i :")
            x_val = Dict((i,j,k) => sol.x[(i,j),k] for (i,j) in A for k in a_nodes)
            y_val = sol.y
            all_cover_set = find_cover_set(A, succ_A, x_val, y_val, a_nodes, d, max_flt, L_M)      
            for element in all_cover_set 
                #println(element)
                cover_set = element[1]
                len_cvs = element[2]
                ac = element[3]
                if cover_set in added_cover_sets
                    #println("DÉJÀ TROUVÉE")
                continue
                else 
                    arc_cover_set = [(cover_set[i],cover_set[i+1]) for i in 1:len_cvs-1]
                    mtn_possibilities = intersect(A_M, arc_cover_set)
                    j_poss = unique([arc[2] for arc in mtn_possibilities])
                    compt_cst += 1
                    @constraint(master_model, sum(x[a,ac] for a in arc_cover_set) <= (len_cvs - 2)*(1 + sum(y[j,ac] for j in j_poss)))
                    push!(added_cover_sets, cover_set)
                end
            end  
        end
        #= x_core = Dict((i,j) => all_solutions[1].x[(i,j)] for (i,j) in A)
        y_core = Dict(j => all_solutions[1].y[j] for j in L_M) =#
        println("NUMBER OF COVER INEQUALITIES : $compt_cst")

        set_optimizer_attribute(master_model, "PoolSearchMode", 0)  # 0 = désactivé
        #set_optimizer_attribute(master_model, "PoolSolutions", 1)   # Retour à la valeur par défaut
    end 
     =#
    prep1_time, prep2_time = 0.0, 0.0
    #prep1_time = 0.0

    #Construction des ensembles utiles pour la suite 
    A_theta2 = [(i,j) for (i,j) in A if j != "t"]
    A_theta3 = [(i,j) for (i,j) in A_M_bar if j != "t"]
    eq_nodes = union(["s"], a_nodes)
    instance["A_theta2"] = A_theta2
    instance["A_theta3"] = A_theta3
    instance["eq_nodes"] = eq_nodes

    #Initialisation des variables du problème complet
    rho_val, lambda_val, phi_val = Dict(j => 0.0 for j in L_M), Dict(j => 0.0 for j in L_M), Dict(j => 0.0 for j in L_M)
    u_val, v_val, w_val = Dict(j => 0.0 for j in V), Dict(j => 0.0 for j in V), Dict(j => 0.0 for j in V)

    #Création du modèle du sous-problème dual
    spdual_model = build_spdual(env, instance)
    
    #Création du modèle du sous-problème dual auxiliaire pour les coupes de Pareto
    aux_spdual, ref_map = copy_model(spdual_model)
    set_optimizer(aux_spdual, Gurobi.Optimizer)
    set_optimizer_attribute(aux_spdual, "OutputFlag", 0)
    var_names = [:theta1, :theta2, :sigma, :theta3, :theta4, :pi1, :alpha1, :alpha2]
    for var_name in var_names
        aux_spdual[var_name] = ref_map[spdual_model[var_name]]
    end
    aux_spdual.ext[:c5_ref] = nothing

    #Les compteurs 
    nbcp_fais, nbcp_opti, nbcp_cover, nb_iterations = 0, 0, 0, 0
    acd_sptime = 0
    last_sp_obj = 0
    compt = 0
    #Le point intérieur
    x_core, y_core = master_LP(env, A, A_T, L_M, V_wt_st, a_nodes, succ_A, pred_A, pred_A_M, acritik_set, ac_critique)
    #x_core , y_core = find_interior_point1(x_core, y_core, A, L_M, a_nodes)
    #feas = MP_is_feasible(x_core, y_core, master_model, a_nodes, A, L_M)
    #= println("\n\n\nPOINT INTÉRIEUR TROUVÉ POUR LE MP ET VÉRIFICATION DE SA FAISABILITÉ : ")
    if !feas
        println("PAS FAISABLE")
    end 
     =#
    #x_core, y_core = nothing, nothing
    mtn_node_inf = Dict((k,j) => 0 for k in acritik_set for j in L_M)
    function my_callback(cb_data)
        status = callback_node_status(cb_data, master_model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end
        nb_iterations += 1

        # Récupérer la solution courante
        x_val = Dict((i,j,k) => round(callback_value(cb_data, x[(i,j),k])) for (i,j) in A for k in a_nodes)
        y_val = Dict((j,k) => round(callback_value(cb_data, y[j,k])) for j in L_M for k in a_nodes)
        Z_val = callback_value(cb_data, Z)
        
        ac_paths, ac_mtn, ac_cumul = build_all_paths(x_val, y_val, d, A, L_M, a_nodes)
        #Résolution du sous-problème dual
        result = solve_spdual(x_val, y_val, instance, spdual_model)
        spd_obj = result.obj
        
        acd_sptime += result.time
        #= out_result = sp_analytic(ac_paths, ac_mtn, d, L_M, max_flt)
        out_srho = out_result.s_rho
        for ac in acritik_set
            if ac_mtn[ac] != nothing && out_srho[ac_mtn[ac]] >= 0
                mtn_node_inf[(ac,ac_mtn[ac])] = 0
            end
        end 
         =#
        if result.status == MOI.OPTIMAL
            #= su_temp, srho_temp, run_time = sp_analytic(ac_paths, ac_mtn, d, L_M, max_flt)
            obj_anal = 0
            for j in L_M 
                obj_anal += srho_temp[j] 
            end  
            println("\nOBJ DU SP_DUAL = ", result.obj, " | TEMPS = ", result.time)
            println("OBJ ANALYTIC = ", obj_anal, " | TEMPS = ", run_time)
            =#

            #= out_result = sp_analytic(ac_paths, ac_mtn, d, L_M, max_flt)
            out_sol = out_result.sol
             =##= println()
            println(out_result.sol)
            =#
            #= println("\nITERATIONS : $nb_iterations")
            println("OBJ DU SP_DUAL = ", result.obj, " | Z_VAL = ", Z_val)
            =#
            gap = spd_obj - Z_val
            #=  
            for k in a_nodes
                for j in L_M
                    y_core[j,k] = 0.7*y_core[j,k] + 0.3*y_val[(j,k)]
                end
            end 
             =#
            #Résolution du sous-problème auxiliaire dual pour la coupe de pareto
            result_aux = solve_aux_spdual(x_val, y_val, x_core, y_core, instance, aux_spdual, spd_obj)
            #= println("OBJ DU SP_DUAL PARETO = ", result_aux.obj, " | Z_VAL = ", Z_val)
            VAL = compare_val(result, x_core, y_core, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
            println("COMPARAISON DES DEUX SOLUTIONS : ", VAL)
             =#
            result = result_aux
            acd_sptime += result.time
        
            if gap > 1e-6
                nbcp_opti += 1
                #= println("\nITERATIONS : $nb_iterations")  
                println("\nSOLUTION FAISABLE") 
                print_x_y(x_val, y_val, A, L_M, d, acritik_set) =#
                 
                #Contruction et ajout de la coupe d'optimalité
                coupe = build_opti_cut(result, x, y, Z, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
                #println(coupe)
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), coupe)
                # Génération de coupes dans le voisinage de y
                #= println("\nITERATIONS : $nb_iterations")  
                println("\nSOLUTION FAISABLE") 
                print_x_y(x_val, y_val, A, L_M, d, acritik_set)
                 =#
                y_list = y_neighborhood(ac_paths, ac_mtn, y_val, acritik_set, A, L_M, A_M)
                for y_sol in y_list
                    result = solve_spdual(x_val, y_sol, instance, spdual_model)
                    acd_sptime += result.time
                    if result.status == MOI.OPTIMAL
                        coupe = build_opti_cut(result, x, y, Z, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
                        MOI.submit(master_model, MOI.LazyConstraint(cb_data), coupe)
                        nbcp_opti += 1
                    elseif result.status == MOI.DUAL_INFEASIBLE
                        coupe = build_fais_cut(result, x, y, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
                        MOI.submit(master_model, MOI.LazyConstraint(cb_data), coupe)
                        nbcp_fais += 1   
                    end    
                end
                #= if compt == 0
                    for k in acritik_set
                        if ac_mtn[k] == nothing
                            cst = @build_constraint(sum(y[j,k] for j in L_M) == 0)
                            #println(cst)
                            MOI.submit(master_model, MOI.LazyConstraint(cb_data), cst)
                            compt += 1
                        end 
                    end
                end =#
            else 
                last_sp_obj = spd_obj
                spd_c1, spd_c2, spd_c3, spd_c4 = spdual_model[:spd_c1], spdual_model[:spd_c2], spdual_model[:spd_c3], spdual_model[:spd_c4]

                for j in L_M
                    rho_val[j] = abs(dual(spd_c1[j]))
                    u_val[j] = abs(dual(spd_c3[j]))
                end
                for j in union(a_nodes, ["s"])
                    u_val[j] = abs(dual(spd_c4[j]))
                end
                for j in setdiff(V, union(a_nodes, L_M, ["s"]))
                    u_val[j] = abs(dual(spd_c2[j]))
                end
            end
            
        elseif result.status == MOI.DUAL_INFEASIBLE 
            #= println("\nITERATIONS : $nb_iterations (SOLUTION INFAISABLE)") 
            print_x_y(x_val, y_val, A, L_M, d, acritik_set)
             =#
            #print_information(result, instance, A_theta2, A_theta3, eq_nodes)
            #println("\tSTATUS DU SPDUAL : $(result.status)")
            #= for ac in acritik_set
                if ac_mtn[ac] != nothing && out_srho[ac_mtn[ac]] < 0
                    println("Aircraft = $ac | mtn_node = $(ac_mtn[ac])")
                    mtn_node_inf[(ac,ac_mtn[ac])] +=1
                    println(mtn_node_inf[(ac,ac_mtn[ac])])
                    if mtn_node_inf[(ac,ac_mtn[ac])] >= 10
                        cst = @build_constraint(y[ac_mtn[ac],ac] == 0)
                        #println(cst)
                        MOI.submit(master_model, MOI.LazyConstraint(cb_data), cst)
                    end
                    #option_after = maintenance_option(ac_mtn[ac], instance)
                    #println(option_after)
                    #= for j in option_after
                        cst = @build_constraint(y[j,ac]  == 0)
                        #println(cst)
                        MOI.submit(master_model, MOI.LazyConstraint(cb_data), cst)
                    end =# 
                end 
            end 
            =#
            nbcp_fais += 1   
            coupe = build_fais_cut(result, x, y, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
            MOI.submit(master_model, MOI.LazyConstraint(cb_data), coupe)
            all_cover_set = find_cover_set(ac_paths, y_val, acritik_set, succ_A, d, max_flt, L_M)  
            #all_cover_set = find_cover_set(A, succ_A, x_val, y_val, a_nodes, d, max_flt, L_M)  
            for (cover_set, len_cvs, ac) in all_cover_set
                #push!(added_cover_sets, cover_set)
                #= println("\nAdding cover set for aircraft $ac : ")
                println(cover_set)
                 =#arc_cover_set = [(cover_set[i],cover_set[i+1]) for i in 1:len_cvs-1]
                mtn_possibilities = intersect(A_M, arc_cover_set)
                j_poss = unique([arc[2] for arc in mtn_possibilities])
                cst1 = @build_constraint(sum(x[a,ac] for a in arc_cover_set) <= (len_cvs - 2)*(1 + sum(y[j,ac] for j in j_poss)))
                
                #cst1 = @build_constraint(sum(sum(x[(i,j),ac] for i in get(pred_A, j, [])) for j in cover_set) <= (len_cvs - 2)*(1 + sum(y[j,ac] for j in intersect(L_M, cover_set))))
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), cst1)
                nbcp_cover += 1
            end   
           
        else
            println("Erreur: Statut du sous-problème = $(result.status)")
        end
        
    end
  
    # Enregistrer le callback
    set_attribute(master_model, MOI.LazyConstraintCallback(), my_callback)
    # Résoudre
    optimize!(master_model)
    println("LAST SP OBJ = ", last_sp_obj)

    status = termination_status(master_model)
    nb_coupes = nbcp_fais + nbcp_opti + nbcp_cover
    x_val = value.(x)
    y_val = value.(y)
    #print_x_y(output_file, x_val, y_val, A, L_M, d)
    obj_val = objective_value(master_model)
    time = round(solve_time(master_model), digits = 2)
    gap = round(relative_gap(master_model)*100, digits = 2)
    nbr_nodes =  MOI.get(master_model, MOI.NodeCount())
    dual_obj = round(objective_bound(master_model); digits = 2)
    time = round(solve_time(master_model), digits = 2)
    acd_sptime = round(acd_sptime, digits = 2)
    mp_time = round(time - acd_sptime, digits = 2)
    nbr_mtn = round(Int, sum(y_val))
    prep1_time = round(prep1_time, digits = 2)
    prep2_time = round(prep2_time, digits = 2)
    final_time  = round(time + prep1_time + prep2_time, digits = 2)
    mtn_stations_used = []
    active_j = Set(j for j in L_M if sum(value(y[j,k]) for k in a_nodes) >= 0.9)
    remaining_stations = setdiff(Set(mtn_stations), mtn_stations_used)

    for ms in remaining_stations
        if !isempty(intersect(Set(L_MS[ms]), active_j))
            push!(mtn_stations_used, ms)
        end
    end
    
    nbr_sts_used = round(Int, length(mtn_stations_used))

    output_results = Dict( "instance" => instance, "obj" => obj_val, "x" => x_val, "y" => y_val, 
                    "u" => u_val, "rho" => rho_val, "v" => v_val, "lambda" => lambda_val,"w" => w_val, 
                    "phi" => phi_val, "gap" => gap, "nbr_nodes" => nbr_nodes, "time" => final_time, "mp_time" => mp_time,
                    "sp_time" => acd_sptime, "dual_obj" => dual_obj, "nbr_mtn" => nbr_mtn,  
                    "nbr_sts_used" => nbr_sts_used, "mtn_stations_used" => mtn_stations_used, 
                    "status" => status, "arc_reduc" => instance["arc_reduc"], 
                    "node_reduc" => instance["node_reduc"], "nbr_iterations" => nb_iterations, 
                    "nbr_cuts" => nb_coupes, "nbr_feasibility_cuts" => nbcp_fais, 
                    "nbr_optimality_cuts" => nbcp_opti, "nbr_cover_cuts" => nbcp_cover, 
                    "prep1_time" => prep1_time, "prep2_time" => prep2_time)
    
    if MTN_CAP
        sz = JuMP.value.(z)
        output_results["z"] = sz
    end

    return output_results
end