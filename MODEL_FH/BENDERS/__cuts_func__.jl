
function build_opti_cut(sub_paths, x, Z, rho_val, fix_dict, added_cuts)
    cut_dict = Dict("opti" => Vector(), "fix" => Vector())
    sum_1 = AffExpr(0.0) 
    for (aircraft, (sub_path, len_sp, mtn_node)) in sub_paths
        if rho_val[mtn_node] > 0
            sum_1 = rho_val[mtn_node]*(1 - sum((1 - x[(sub_path[i],sub_path[i+1])]) for i in 1:len_sp-1))
            cut = @build_constraint(sum_1 <= Z[aircraft]) 
            push!(cut_dict["opti"], cut)
            #push!(added_cuts, sub_path)
        elseif rho_val[mtn_node] == 0
            fix_dict[aircraft] += 1
            if fix_dict[aircraft] == 5
                cut = @build_constraint(sum((1 - x[(sub_path[i],sub_path[i+1])]) for i in 1:len_sp-1) <= 0)
                push!(cut_dict["fix"], cut)
            end 
        end  
    end
    return cut_dict
end 

function build_opti_cut_agg(sub_paths, x, Z, obj_sp, rho_val)
    summ = AffExpr(0.0) 
    cut_list = []
    for (aircraft, (sub_path, len_sp, mtn_node)) in sub_paths
        if len_sp != 0
            summ += sum((1 - x[(sub_path[i],sub_path[i+1])]) for i in 1:len_sp-1)
        end
        if rho_val[mtn_node] == 0
            #println(aircraft)
            cut = @build_constraint(sum((1 - x[(sub_path[i],sub_path[i+1])]) for i in 1:len_sp-1) <= 0)
            #cut = @build_constraint(Z[aircraft] <= 0)
            push!(cut_list, cut) 
        end
    end
    cut = @build_constraint(obj_sp*(1-summ) <= Z) 
    push!(cut_list, cut) 

    return cut_list 
end 


function build_feas_cut(infeas_subpaths, x, d, A_M)
    cut_list = []
    for (_, (sub_path, len_sp)) in infeas_subpaths
        temp = sum(d[i] for i in sub_path[1:end-1])
        if temp == 6000 && (sub_path[end-1], sub_path[end]) in A_M
            println("\nFeas cuts : ", temp)
            println(sub_path)
        end 
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


function format_gap(UB, LB)
    if UB == Inf || UB == -Inf
        return "N/A"
    end
    if UB == LB  # inclut le cas 0/0
        return "0.0%"
    end
    gap = (UB - LB)/UB*100
    return "$(round(gap, digits=2))%"
end


