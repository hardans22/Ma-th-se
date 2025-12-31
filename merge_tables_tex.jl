function simple_merge_tables(file1, file2, output_file, new_caption)
    """Version simplifiée pour tableaux basiques"""
    
    lines1 = readlines(file1)
    lines2 = readlines(file2)
    
    merged = String[]
    
    for (i, (l1, l2)) in enumerate(zip(lines1, lines2))
        if occursin("\\begin{tabular}", l1)
            # Extraire les formats de colonnes
            format1 = match(r"\\begin\{tabular\}\{([^}]+)\}", l1).captures[1]
            format2 = match(r"\\begin\{tabular\}\{([^}]+)\}", l2).captures[1]
            # Combiner proprement
            push!(merged, "\\begin{tabular}{$(format1)$(format2[2:end])}")
        elseif occursin("\\end{tabular}", l1)
            push!(merged, l1)
        elseif occursin("\\\\hline", l1)
            clean_l1 = replace(l1, r"\\\\\\hline?\s*$" => "")
            clean_l2 = lstrip(replace(l2, r"^[^&]*" => ""))
            push!(merged, clean_l1 * clean_l2)
        elseif occursin("\\caption", l1)
            push!(merged, "\\caption{$(new_caption)}")
        else
            # Combiner les lignes de données
            clean_l1 = replace(l1, r"\\\\?\s*$" => "")
            clean_l2 = lstrip(replace(l2, r"^[^&]*" => ""))
            push!(merged, clean_l1 * clean_l2)
        end
    end
    
    open(output_file, "w") do f
        for line in merged
            println(f, line)
        end
    end
end

# Utilisation
file1 = "RESULTS/A_MTN_2_SMALL_INSTANCES_/summary_FC_MAX_DY_MIN.txt"
file2 = "RESULTS/A_MTN_2_SMALL_INSTANCES_/summary_FC_MAX_DY_MAX.txt"
file3 = "RESULTS/A_MTN_2_SMALL_INSTANCES_/summary_FC_MAX_DY_MIN_MAX.txt"
new_caption = "Résumé des résultats - summary\\_FC\\_MAX\\_DY\\_MIN\\_MAX"


simple_merge_tables(file1, file2, file3, new_caption)