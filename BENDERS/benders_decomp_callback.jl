include("functions_for_benders.jl")

function benders_decomp(env, instance_file, output_file, nbr_thread, silent, time_limit)
    instance = build_graph(instance_file, false)
    nbr_K = instance["number_of_aircrafts"]
    max_flt = instance["maximum_flying_time"]
    d = instance["d"]
    A_S, A_T = instance["A_S"], instance["A_T"]
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
    @variable(master_model, x[k in A], Bin)
    @variable(master_model, y[j in L_M, k in a_nodes], Bin)
    @variable(master_model, Z >= 0)

    A_M_bar_without_t = [a for a in A_M_bar if a[2] != "t"]

    # Contraintes
    @constraint(master_model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(master_model, sum(x[i] for i in A_T) == nbr_K)
    @constraint(master_model, mp_c1[i in fl_nodes], sum(x[(i,j)] for j in get(succ_A, i, [])) == 1)
    @constraint(master_model, mp_c2[i in V_wt_st], sum(x[(j,i)] for j in get(pred_A, i, [])) == sum(x[(i,j)] for j in get(succ_A, i, [])))
    @constraint(master_model, mp_c3[j in L_M], sum(y[j,k] for k in a_nodes) <= sum(x[(i, j)] for i in get(pred_A_M, j, [])))
    
    #@constraint(master_model, sum(y[j,k] for j in L_M, k in a_nodes) <= ac_critique)
    @constraint(master_model, mp_c4[k in acritik_set], sum(y[j,k] for j in L_M) <= 1)
    @constraint(master_model, mp_c5[j in L_M, k in setdiff(a_nodes, acritik_set)], sum(y[j,k] for j in L_M) == 0)
    #write_to_file(master_model, "master_model.lp")
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
            ac = element[3]
            if cover_set in added_cover_sets
                #println("DÉJÀ TROUVÉE")
                continue
            else 
                arc_cover_set = [(cover_set[i],cover_set[i+1]) for i in 1:len_cvs-1]
                mtn_possibilities = intersect(A_M, arc_cover_set)
                j_poss = unique([arc[2] for arc in mtn_possibilities])
                @constraint(master_model, sum(x[a] for a in arc_cover_set) <= (len_cvs - 2)*(1 + sum(y[j,ac] for j in j_poss)))
                push!(added_cover_sets, cover_set)
            end
        end
    end
 =#
    prep2_time = @elapsed begin
        set_optimizer_attribute(master_model, "PoolSearchMode", 2)  # Mode exhaustif
        set_optimizer_attribute(master_model, "PoolSolutions",10)  # Nombre max de solutions à stocker
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
            x_val = sol.x
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
                    @constraint(master_model, sum(x[a] for a in arc_cover_set) <= (len_cvs - 2)*(1 + sum(y[j,ac] for j in j_poss)))
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
  
    prep1_time = 0.0
    
    A_theta2 = [(i,j) for (i,j) in A if j != "t"]
    A_theta3 = [(i,j) for (i,j) in A_M_bar if j != "t"]
    eq_nodes = union(["s"], a_nodes)
    instance["A_theta2"] = A_theta2
    instance["A_theta3"] = A_theta3
    instance["eq_nodes"] = eq_nodes
    
    spdual_model = build_spdual(env, instance)
    sp_model = build_sp(env, instance)
    
    rho_val, lambda_val, phi_val = Dict(j => 0.0 for j in L_M), Dict(j => 0.0 for j in L_M), Dict(j => 0.0 for j in L_M)
    u_val, v_val, w_val = Dict(j => 0.0 for j in V), Dict(j => 0.0 for j in V), Dict(j => 0.0 for j in V)

    aux_spdual, ref_map = copy_model(spdual_model)
    set_optimizer(aux_spdual, Gurobi.Optimizer)
    set_optimizer_attribute(aux_spdual, "OutputFlag", 0)
    # Liste des noms de variables à mapper
    var_names = [:theta1, :theta2, :sigma, :theta3, :theta4, :pi1, :alpha1, :alpha2]

    # Mapper automatiquement
    for var_name in var_names
        aux_spdual[var_name] = ref_map[spdual_model[var_name]]
    end
    aux_spdual.ext[:c5_ref] = nothing

    nbcp_fais, nbcp_opti, nbcp_cover, nb_iterations = 0, 0, 0, 0
    acd_sptime = 0
    last_sp_obj = 0

    added_cover_sets = Set{Vector{String}}()  # ou le type de vos nœuds
    node_count = 0

    x_core = nothing
    y_core = nothing
    
    max_history = 5

    function my_callback(cb_data)
        status = callback_node_status(cb_data, master_model)
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end
        # Obtenir le numéro de nœud via l'API Gurobi directement
        #= node_count = Ref{Cdouble}()
        Gurobi.GRBcbget(cb_data, GRB_CB_MIPSOL, GRB_CB_MIPSOL_NODCNT, node_count)
        
        println("Nœud actuel : ", node_count[])
         =#
        nb_iterations += 1

        # Récupérer la solution courante
        x_val = Dict((i,j) => callback_value(cb_data, x[(i,j)]) for (i,j) in A)
        y_val = Dict((j,k) => callback_value(cb_data, y[j,k]) for j in L_M for k in a_nodes)
        Z_val = callback_value(cb_data, Z)
        
        result = solve_spdual(x_val, y_val, instance, spdual_model)
        #= result_sp = solve_sp(sp_model, x_val, y_val)
        if result.status == MOI.DUAL_INFEASIBLE && result_sp.status != MOI.INFEASIBLE
            println("PROBLEMMMMMMMMME")
        end =#
        spd_obj = result.obj
        #= n_solutions = result_count(spdual_model)
        println("Nombre de solutions trouvées : $n_solutions")
        =#
        acd_sptime += result.time
        #println(result.status)
        if result.status == MOI.OPTIMAL
            #= println("\nITERATIONS : $nb_iterations")
            print_x_y(x_val, y_val, A, L_M, d, a_nodes)
            println("OBJ DU SP_DUAL = ", result.obj, " | Z_VAL = ", Z_val)
          =#
            #= coupe = build_opti_cut(result, x, y, Z, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
            println(coupe)
             =##print_information(result, instance, A_theta2, A_theta3, eq_nodes)
            
            gap = spd_obj - Z_val
            if x_core == nothing 
                x_core = copy(x_val)
                y_core = copy(y_val)
            end 
            #= VAL = compare_val(result, x_core, y_core, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
            println("VAL NORMAL : $VAL")
             =#
            #= for (i,j) in A 
                x_core[(i,j)] = 0.5*x_core[(i,j)] + 0.5*x_val[(i,j)]
            end 
            for j in L_M
                y_core[j] = 0.5*y_core[j] + 0.5*y_val[j]
            end
             =#result_aux = solve_aux_spdual(x_val, y_val, x_core, y_core, instance, aux_spdual, spd_obj)
            #println("OBJ DU SP_DUAL PARETO = ", result_aux.obj, " | Z_VAL = ", Z_val)
            #= cut_comparison(result, result_aux, x_val, y_val, instance, A_theta2, A_theta3)
            coupe = build_opti_cut(result_aux, x, y, Z, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
            println(coupe)
             =#
            result = result_aux
            #print_information(result, instance, A_theta2, A_theta3, eq_nodes)
            acd_sptime += result.time
            #println("\tSTATUS DU SPDUAL : $(result.status)")
            if gap > 1e-6
                nbcp_opti += 1
                
                #result = cut_comparison(result, result_aux, x_val, y_val, instance, A_theta2, A_theta3)
                coupe = build_opti_cut(result, x, y, Z, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), coupe)
                y_list = y_neighborhood(x_val, y_val, acritik_set, A, L_M, A_M)
                #println("NOMBRE DE SOLUTIONS DANS LE VOISINAGE DE Y : ", length(y_list))
                 for y_sol in y_list
                    result = solve_spdual(x_val, y_sol, instance, spdual_model)
                    acd_sptime += result.time
                    if result.status == MOI.DUAL_INFEASIBLE
                        #nbcp_fais += 1   
                        coupe = build_fais_cut(result, x, y, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
                    elseif result.status == MOI.OPTIMAL
                        #nbcp_opti += 1 
                        coupe = build_opti_cut(result, x, y, Z, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
                    end 
                    MOI.submit(master_model, MOI.LazyConstraint(cb_data), coupe)
                end 
            else 
                #println("\tELSE")
                last_sp_obj = spd_obj
                #println("OBJ DU SP_DUAL = ", result.obj, " | Z_VAL = ", Z_val, " | GAP = ", gap)
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
            #println("\tSTATUS DU SPDUAL : $(result.status)")
            print_x_y(x_val, y_val, A, L_M, d, a_nodes)

            nbcp_fais += 1   
            coupe = build_fais_cut(result, x, y, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
            MOI.submit(master_model, MOI.LazyConstraint(cb_data), coupe)
            #=x_neighbors = x_neighborhood(x_val, fl_nodes, pred_A, succ_A, V)
            #println(length(x_neighbors))
             for x_new in x_neighbors
                result = solve_spdual(x_new, y_val, instance, spdual_model)
                if result.status == MOI.DUAL_INFEASIBLE
                    #nbcp_fais += 1   
                    acd_sptime += result.time
                    coupe = build_fais_cut(result, x, y, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
                    MOI.submit(master_model, MOI.LazyConstraint(cb_data), coupe)
                else 
                    acd_sptime += result.time
                    coupe = build_opti_cut(result, x, y, Z, max_flt, d, f, A_M, A_theta2, A_theta3, L_M, a_nodes, V)
                    MOI.submit(master_model, MOI.LazyConstraint(cb_data), coupe)
                end
            end  
             =##println("\tCoupe ajoutée : $(coupe)")
            #print_x_y(x_values, y_values, A, L_M, d)
            
            all_cover_set = find_cover_set(A, succ_A, x_val, y_val, acritik_set, d, max_flt, L_M)  
            #best_mtn_set = find_mtn_set(all_cover_set, A_M, d, max_flt)
            #println("\nBEST MTN_SET = ", best_mtn_set.bnd_val)
            
            #println("NOMBRE DE COVER SETS : ", length(all_cover_set))          
            for (cover_set, len_cvs, ac) in all_cover_set
                #println(cover_set)
                
                if cover_set in added_cover_sets
                    println("\n\n\nCOVER DÉJÀ TROUVÉE PLUS HAUT \n\n")
                    continue
                end 
                push!(added_cover_sets, cover_set)
                arc_cover_set = [(cover_set[i],cover_set[i+1]) for i in 1:len_cvs-1]
                mtn_possibilities = intersect(A_M, arc_cover_set)
                j_poss = unique([arc[2] for arc in mtn_possibilities])
                cst1 = @build_constraint(sum(x[a] for a in arc_cover_set) <= (len_cvs - 2)*(1 + sum(y[j,ac] for j in j_poss)))
                #println(cst1)
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), cst1)
                nbcp_cover += 1
            end  
           
        else
            println("Erreur: Statut du sous-problème = $(result.status)")
        end
        
    end
    #write_to_file(master_model, "mp.lp")
    # Enregistrer le callback
    set_attribute(master_model, MOI.LazyConstraintCallback(), my_callback)
    #set_attribute(master_model, MOI.RawOptimizerAttribute("PreCrush"), 1)
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