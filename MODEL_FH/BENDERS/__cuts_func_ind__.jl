function build_opti_cut_ind(aircraft, irr_path, x, Z, result)
    cut = nothing
    sum_1 = AffExpr(0.0)
   #=  if aircraft =="N122"
        println("Affiché dans la fonction de coupe :", result.obj)
        #= temp = sum([d[x] for x in aircraft_paths[ac].path])
        println(temp)
         =#
    end =#
    if irr_path !== nothing
        sub_path, len_sp = irr_path[1], irr_path[2]
    else
        return nothing
    end
    #println(irr_path)
    sum_1 = result.obj*(1 - sum((1 - x[(sub_path[i],sub_path[i+1])]) for i in 1:len_sp-1))
    cut = @build_constraint(sum_1 <= Z[aircraft])
    return cut
end 


function build_fix_cut(aircraft, ac_path, irr_path, x, fix_dict)
    cut = nothing
    if irr_path !== nothing
        sub_path, len_sp = irr_path[1], irr_path[2]
    else
        sub_path, len_sp = ac_path.path, ac_path.len_path 
    end
    fix_dict[aircraft] += 1
    if fix_dict[aircraft] == 1
        cut = @build_constraint(sum((1 - x[(sub_path[i],sub_path[i+1])]) for i in 1:len_sp-1) <= 0)
    end 
    return cut
end 


function build_feas_cut_ind(irr_path, x)
    cut = nothing
    #println(irr_path           )
    sub_path, len_sp = irr_path[1], irr_path[2]
    cut  = @build_constraint(sum(x[(sub_path[i],sub_path[i+1])] for  i in 1:len_sp-1) <= len_sp - 2)    
    return cut                                 
end       

