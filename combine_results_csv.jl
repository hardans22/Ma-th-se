using CSV, DataFrames, Statistics, Glob

FC = true 
DY = true

#= function combine_csv_FH(folder, ind)
    for csv_type in ["summary", "aircraft"]
        files = filter(f -> endswith(f, ".csv") && occursin(csv_type, lowercase(f)), 
                      readdir(folder, join=true))
        df_combined = vcat([CSV.read(f, DataFrame) for f in files]...)
        
        if ind == "FC" && csv_type == "summary"
            select!(df_combined, Not(:fh_day))
        elseif ind == "DY" && csv_type == "summary"
            select!(df_combined, Not(:fh_tk))
        end
        
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
                    df_combined[!, col] = round.(df_combined[!, col], digits=1)
                end
            end
        end
        
        name = split(folder, "/")[end-1]
        path = dirname(dirname(folder)) * "/" * csv_type * "_" * name * ".csv"
        
        CSV.write(path, df_combined)
        println("Fichier créé : $path")
    end
end

if FC  
    folder_1 = "RESULTS/A_MTN_2_SMALL_INSTANCES_/FH_FC_"
    for type in ["MIN", "MAX", "MEDIAN", "MEAN"]
        folder = folder_1 * type * "/"
        combine_csv(folder, "FC")
    end    
end

if DY 
    folder_1 = "RESULTS/A_MTN_2_SMALL_INSTANCES_/FH_DY_"
    for type in ["MIN", "MAX", "MEDIAN", "MEAN"]
        folder = folder_1 * type * "/"
        combine_csv(folder, "DY")
    end
end =#

function combine_csv_FC(folder)
    for csv_type in ["summary", "aircraft"]
        files = filter(f -> endswith(f, ".csv") && occursin(csv_type, lowercase(f)), 
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
                    df_combined[!, col] = round.(df_combined[!, col], digits=1)
                end
            end
        end
        
        name = split(folder, "/")[end]
        path = dirname(folder) * "/" * csv_type * "_" * name * ".csv"
        
        CSV.write(path, df_combined)
        println("Fichier créé : $path")
    end
end


folder_1 = "RESULTS/A_MTN_2_SMALL_INSTANCES_/"

fc_folders = glob("FC*", folder_1)

for folder in fc_folders
    combine_csv_FC(folder)
end 