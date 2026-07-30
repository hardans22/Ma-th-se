using CSV, DataFrames, Statistics, Glob

function combine_csv(folder, type)
    for csv_type in ["summary"]
        files = filter(f -> endswith(f, type) && occursin(csv_type, lowercase(f)), 
                      readdir(folder, join=true))
        df_combined = vcat([CSV.read(f, DataFrame) for f in files]...)
        
        function split_code(s::AbstractString)
            parts = split(s, "_")
            prefix = join(parts[1:end-1], "_")
            number = parse(Int, parts[end])
            return prefix, number
        end
        
        sort!(df_combined, :Instances, by = s -> split_code(s))
        
        # Ajouter les moyennes par type d'instance (pour summary uniquement)
        if csv_type == "summary"
            # Extraire le préfixe de chaque instance
            df_combined.instance_type = [split_code(inst)[1] for inst in df_combined.Instances]
            
            # Grouper par type d'instance et calculer les moyennes
            grouped = groupby(df_combined, :instance_type)
            
            df_with_means = DataFrame()
            for group in grouped
                # Ajouter les lignes du groupe
                append!(df_with_means, select(group, Not(:instance_type)), promote=true)
                
                # Créer la ligne de moyenne
                mean_row = DataFrame()
                mean_row.Instances = ["AVERAGE"]
                
                # Calculer la moyenne pour chaque colonne numérique
                for col in names(group)
                    if col != "Instances" && col != "instance_type"
                        if eltype(group[!, col]) <: Number
                            mean_row[!, col] = [mean(skipmissing(group[!, col]))]
                        else
                            mean_row[!, col] = [""]
                        end
                    end
                end
                
                append!(df_with_means, mean_row, promote=true)
            end
            
            df_combined = df_with_means
            
            # Arrondir les colonnes décimales à 2 chiffres
            for col in names(df_combined)
                if col != "Instances" && eltype(df_combined[!, col]) <: AbstractFloat
                    digits = occursin(r"time|temps|cpu"i, col) ? 4 : 1
                    df_combined[!, col] = round.(df_combined[!, col], digits=digits)
                end
            end
        end
        
        name = split(folder, "/")[end]
        path = dirname(folder) * "/" * csv_type * "_" * name * type
        
        CSV.write(path, df_combined)
        println("Fichier créé : $path")
    end
end


for fold in ["A_MTN_3/"]
    println("FOLD = ", fold)
    folder = "/home/dansou/Téléchargements/RESULTS_MILP/"*fold
    type = "MILP.csv"
    combine_csv(folder,type)
end