function build_opti_cut_ind(aircraft, irr_path, x, Z, result, fix_dict, added_cuts)
    cut = nothing
    cut_type = nothing 
    sum_1 = AffExpr(0.0)
    sub_path, len_sp = irr_path[1], irr_path[2]
    if result.obj > 0
        sum_1 = result.obj(1 - sum((1 - x[(sub_path[i],sub_path[i+1])]) for i in 1:len_sp-1))
        cut = @build_constraint(sum_1 <= Z[aircraft]) 
        cut_type = "opti"
    elseif result.obj == 0
        fix_dict[aircraft] += 1
        if fix_dict[aircraft] == 5
            cut = @build_constraint(sum((1 - x[(sub_path[i],sub_path[i+1])]) for i in 1:len_sp-1) <= 0)
            cut_type = "fix"
        end 
    end
    return cut, cut_type
end 

function build_feas_cut_ind(irr_path, x)
    cut = nothing
    sub_path, len_sp = irr_path[1], irr_path[2]
    cut  = @build_constraint(sum(x[(sub_path[i],sub_path[i+1])] for  i in 1:len_sp-1) <= len_sp - 2)    
    return cut
end 

