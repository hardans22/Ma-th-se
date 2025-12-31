using PrettyTables, CSV, DataFrames

function csv_to_latex_pretty(csv_file::String; 
                             caption::String="",
                             label::String="",
                             position::String="htbp",
                             font_size::String="\\small",
                             style::Symbol=:booktabs,
                             format_decimals::Int=2,
                             highlight_mean::Bool=true)
    
    df = CSV.read(csv_file, DataFrame)
    #select!(df, Not(:LB, :Nodes, :Opt, :Nbr_mtn, :Nbr_sts_used))
    col_order = ["Instances", "fh_tk", "fh_day", "Rho", "Lambda", "Phi", "UB", "Feasible", "Gap", "Time"]
    # Garder seulement les colonnes qui existent
    col_order = [col for col in col_order if col in names(df)]
    
    df = select(df, vcat(col_order))
    
    # Renommer les colonnes
    rename!(df, "Rho" => "FH", "Lambda" => "FC", "Phi" => "DY", "Feasible" => "Fsb")
    
    # Récupérer le nom sans extension
    nom = replace(basename(csv_file), ".csv" => "")
    output_file = replace(csv_file, ".csv" => ".txt")
    
    # Caption et label par défaut
    if caption == ""
        caption = "Résultats pour $nom"
    end
    if label == ""
        label = "tab:" * lowercase(replace(nom, "_" => "-"))
    end
    
    # Choisir le style de table
    #= table_style = style == :booktabs ? tf_latex_booktabs : 
                  style == :matrix ? tf_latex_matrix :
                  tf_latex_default =#
    table_style = LatexTableFormat(top_line = "\\hline", header_line = "\\hline",
        mid_line = "\\hline", bottom_line = "\\hline", left_vline = "|", 
        mid_vline = "", right_vline = "|", header_envs = ["textbf"]
    )
    
    # Formater les colonnes numériques
    formatters = ft_printf("%.$(format_decimals)f", 2:ncol(df))
    
    # Highlighter pour les lignes de moyenne
    hl_mean = LatexHighlighter(
        (data, i, j) -> occursin("AVERAGE", string(data[i, 1])),
        ["textbf"])
    
    # Écrire dans un fichier
    open(output_file, "w") do io
        # Ajouter l'environnement table
        println(io, "\\begin{table}[$position]")
        println(io, "\\centering")
        println(io, font_size)
        # Ajouter caption et label
        println(io, "\\caption{$caption}")
        println(io, "\\label{$label}")
        # Générer le tableau avec PrettyTables
        pretty_table(io, df,
                     backend = Val(:latex),
                     tf = table_style,
                     formatters = formatters,
                     highlighters = highlight_mean ? hl_mean : nothing,
                     show_subheader = false,  # Utiliser show_subheader au lieu de nosubheader
                     header = names(df),
                     alignment = :r,  # Colonnes centrées
                     wrap_table = false)  # Ne pas wrapper dans table (on le fait manuellement)
        
       
        println(io, "\\end{table}")
    end
    
    println("✓ Fichier LaTeX créé : $output_file")
end

# Fonction pour créer tous les tableaux
function create_all_latex_tables(base_folder::String; kwargs...)
    csv_files = filter(f -> endswith(f, ".csv") && occursin("summary", lowercase(f)), 
                      readdir(base_folder, join=true))
    
    for csv_file in csv_files
        nom = replace(basename(csv_file), ".csv" => "")
        
        # Personnaliser selon le type de fichier
        if occursin("summary", nom)
            caption = "Résumé des résultats - $nom"
        elseif occursin("aircraft", nom)
            caption = "Détails par aéronef - $nom"
        else
            caption = "Résultats - $nom"
        end
        
        csv_to_latex_pretty(csv_file, caption=caption; kwargs...)
    end
    
    println("\n🎉 Tous les tableaux LaTeX ont été créés!")
end




# Exemple 2 : Créer tous les tableaux
create_all_latex_tables("RESULTS/A_MTN_2_SMALL_INSTANCES_/",font_size="\\small",style=:booktabs,format_decimals=1)


#= files = filter(f -> endswith(f, ".csv") && occursin("summary", lowercase(f)), 
                      readdir(fold_path, join=true))
for file in files
    df = CSV.read(file, DataFrame)
    select!(df, Not(:LB, :Nodes, :Opt, :Nbr_mtn, :Nbr_sts_used))

    # Afficher directement en LaTeX
    println("Fichier traité : $file")

    name = replace(basename(file), ".csv" => ".txt")
    # Ou sauvegarder dans un fichier
    open(fold_path*name, "w") do io
        pretty_table(io, df, backend=Val(:latex),
                    tf=tf_latex_booktabs,  # Style booktabs (plus joli)
                    formatters=ft_printf("%.1f", 2:ncol(df)))  # Format décimal
    end
end =#