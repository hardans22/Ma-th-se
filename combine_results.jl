using DataFrames, XLSX, Glob

function combine_all_results(init_path, output_without, output_with)
    # Initialiser les DataFrames pour stocker tous les résultats
    all_without_data = DataFrame[]
    all_with_data = DataFrame[]
    
    # Créer les listes d'identification
    mw_labels = []
    ms_labels = []
    
    for first_fold in ["A_MTN_2/"]
        current_path = init_path * first_fold
        # Extraire les labels pour identification
        mw_label = replace(first_fold, "/" => "")  # Enlever le "/"
        
        # Traiter les fichiers "without"
        without_files = sort(glob("*_without_*.xlsx", current_path))
        for file in without_files
            println("Processing without file: ", file)
            df = DataFrame(XLSX.readtable(file, "Sheet1"))
            
            # Ajouter les colonnes d'identification au début
            df_with_labels = DataFrame()
            df_with_labels.AC = fill(mw_label, nrow(df))
            
            # Ajouter toutes les autres colonnes
            for col in names(df)
                df_with_labels[!, col] = df[!, col]
            end
            
            push!(all_without_data, df_with_labels)
        end
        
        # Traiter les fichiers "with"
        with_files = sort(glob("*_with_*.xlsx", current_path))
        for file in with_files
            println("Processing with file: ", file)
            df = DataFrame(XLSX.readtable(file, "Sheet1"))
            
            # Ajouter les colonnes d'identification au début
            df_with_labels = DataFrame()
            df_with_labels.AC = fill(mw_label, nrow(df))
            
            # Ajouter toutes les autres colonnes
            for col in names(df)
                df_with_labels[!, col] = df[!, col]
            end
            
            push!(all_with_data, df_with_labels)
        end
    end
    
    # Combiner tous les DataFrames
    if !isempty(all_without_data)
        combined_without = vcat(all_without_data..., cols=:union)
        XLSX.writetable(output_without, combined_without, overwrite=true)
        println("Fichier 'without' créé: ", output_without)
        println("Nombre total de lignes 'without': ", nrow(combined_without))
    end
    
    if !isempty(all_with_data)
        combined_with = vcat(all_with_data..., cols=:union)
        XLSX.writetable(output_with, combined_with, overwrite=true)
        println("Fichier 'with' créé: ", output_with)
        println("Nombre total de lignes 'with': ", nrow(combined_with))
    end
end

# Utilisation
init_path = "/home/harcenage/Téléchargements/RESULTS/"
output_without = init_path * "FINAL_all_results_without_prep.xlsx"
output_with = init_path * "FINAL_all_results_with_prep.xlsx"

combine_all_results(init_path, output_without, output_with)