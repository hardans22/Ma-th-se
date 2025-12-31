using CSV, DataFrames, Dates, JSON, XLSX, Statistics

function process_files(folder_path, file_type)
    # Obtenir tous les fichiers .xlsx
    files = filter(f -> endswith(f, file_type), readdir(folder_path))
    return files
end

function computation(df_flights)
    fh_tk_min = minimum(df_flights.AIR_TIME)
    fh_tk_max = maximum(df_flights.AIR_TIME)
    fh_tk_md = round(Int, median(df_flights.AIR_TIME))
    fh_tk_mn = round(Int, mean(df_flights.AIR_TIME))

    aircraft_day_df = combine(groupby(df_flights, [:TAIL_NUMBER, :DAY]),
        :AIR_TIME => sum => :FLYING_TIME,
        :AIR_TIME => length => :TAKEOFF
    )

    sort!(aircraft_day_df, [:TAIL_NUMBER, :DAY])

    fh_day_min = minimum(aircraft_day_df.FLYING_TIME)
    fh_day_max = maximum(aircraft_day_df.FLYING_TIME)
    fh_day_md = round(Int, median(aircraft_day_df.FLYING_TIME))
    fh_day_mn = round(Int, mean(aircraft_day_df.FLYING_TIME))

    tk_day_min = minimum(aircraft_day_df.TAKEOFF)
    tk_day_max = maximum(aircraft_day_df.TAKEOFF)
    tk_day_md = round(Int, median(aircraft_day_df.TAKEOFF))
    tk_day_mn = round(Int, mean(aircraft_day_df.TAKEOFF))
    
    return (fh_tk_min = fh_tk_min, fh_tk_max = fh_tk_max, fh_tk_md = fh_tk_md, fh_tk_mn = fh_tk_mn,
            fh_day_min = fh_day_min, fh_day_max = fh_day_max, fh_day_md = fh_day_md, fh_day_mn = fh_day_mn,
            tk_day_min = tk_day_min, tk_day_max = tk_day_max, tk_day_md = tk_day_md, tk_day_mn = tk_day_mn)
end

function computation2(df_flights, nbr_ac)
    aircraft_day_df = combine(groupby(df_flights, [:TAIL_NUMBER, :DAY]),
        :AIR_TIME => sum => :FLYING_TIME,
        :AIR_TIME => length => :TAKEOFF
    )

    sort!(aircraft_day_df, [:TAIL_NUMBER, :DAY])
    fh_day = round(Int, sum(aircraft_day_df.FLYING_TIME)/nbr_ac/7)
    tk_day = round(Int, sum(aircraft_day.TAKEOFF)/nbr_ac/7)
    fh_tk = round(sum(Int, aircraft_day.FLYING_TIME)/sum(aircraft_day.TAKEOFF))

    return (fh_day = fh_day, tk_day = tk_day, fh_tk = fh_tk)

end 


#= for ac_critique in [2]  
    folder_path_xlsx ="INSTANCES/instances_xlsx/A_MTN_"*string(ac_critique)*"/"  # ou le chemin vers votre dossier
    folder_path_json ="INSTANCES/instances_json/A_MTN_"*string(ac_critique)*"/"  # ou le chemin vers votre dossier
    xlsx_files = process_files(folder_path_xlsx, "xlsx")
    for file in xlsx_files
        inst_name = split(file, ".")[1]
        file_xlsx = folder_path_xlsx*file
        println("Traitement de: $file_xlsx")
        df_flights = DataFrame(XLSX.readtable(file_xlsx, "Data"))
        
        ind = computation(df_flights)
        coeff_df = DataFrame(FH_TK_MIN = ind.fh_tk_min, FH_TK_MAX = ind.fh_tk_max, FH_TK_MEDIAN = ind.fh_tk_md, FH_TK_MEAN = ind.fh_tk_mn, 
                             FH_DAY_MIN = ind.fh_day_min, FH_DAY_MAX = ind.fh_day_max, FH_DAY_MEDIAN = ind.fh_day_md, FH_DAY_MEAN = ind.fh_day_mn,
                             TK_DAY_MIN = ind.tk_day_min, TK_DAY_MAX = ind.tk_day_max, TK_DAY_MEDIAN = ind.tk_day_md, TK_DAY_MEAN = ind.tk_day_mn,
                             )

        XLSX.openxlsx(file_xlsx, mode="rw") do xf
            # Ajouter une nouvelle feuille
            new_sheet = XLSX.addsheet!(xf, "Coefficients")
            
            # Écrire le DataFrame dans la nouvelle feuille
            XLSX.writetable!(new_sheet, Tables.columntable(coeff_df), write_columnnames=true)
        end
        
        #= # Lire toutes les feuilles existantes
        xf = XLSX.readxlsx(file_xlsx)
        noms_feuilles = XLSX.sheetnames(xf)
        # Préparer toutes les feuilles
        feuilles = Dict()
        for nom in noms_feuilles
            if nom == "Parameters"
                feuilles[Symbol(nom)] = (collect(eachcol(df_param)), names(df_param))
            else
                feuilles[Symbol(nom)] = XLSX.readtable(file_xlsx, nom)
            end
        end
       # Ouvrir le fichier en mode lecture/écriture
        XLSX.openxlsx(file_xlsx, mode="rw") do xf
            sheet = xf["Parameters"]
            
            # Effacer le contenu existant de la feuille
            XLSX.rename!(sheet, "Parameters")
            
            # Écrire le DataFrame
            XLSX.writetable!(sheet, collect(eachcol(df_param)), names(df_param))
        end
         =#

        file_json = folder_path_json*inst_name*".json"
        println("Traitement de: $file_json")
        data = JSON.parsefile(file_json)
        data_dict = Dict(data)
        
        data_dict["fh_tk_min"] = ind.fh_tk_min
        data_dict["fh_tk_max"] = ind.fh_tk_max
        data_dict["fh_tk_median"] = ind.fh_tk_md
        data_dict["fh_tk_mean"] = ind.fh_tk_mn
        data_dict["fh_day_min"] = ind.fh_day_min
        data_dict["fh_day_max"] = ind.fh_day_max
        data_dict["fh_day_median"] = ind.fh_day_md
        data_dict["fh_day_mean"] = ind.fh_day_mn
        data_dict["tk_day_min"] = ind.tk_day_min
        data_dict["tk_day_max"] = ind.tk_day_max
        data_dict["tk_day_median"] = ind.tk_day_md
        data_dict["tk_day_mean"] = ind.tk_day_mn

        # 4. Sauvegarder
        open(file_json, "w") do f
            JSON.print(f, data_dict;)
        end

    end  
 
end
 =#
for ac_critique in [3,5,8]  
    folder_path_xlsx ="INSTANCES/instances_WCTR_xlsx/A_MTN_"*string(ac_critique)*"/"  # ou le chemin vers votre dossier
    folder_path_json ="INSTANCES/instances_WCTR_json/A_MTN_"*string(ac_critique)*"/"  # ou le chemin vers votre dossier
    xlsx_files = process_files(folder_path_xlsx, "xlsx")
    for file in xlsx_files
        inst_name = split(file, ".")[1]
        file_xlsx = folder_path_xlsx*file
        println("Traitement de: $file_xlsx")
        df_ac = DataFrame(XLSX.readtable(file_xlsx, "Aircrafts"))
        
        for row in eachrow(df_ac)
            if row.INIT_FLYING_TIME == 0
                row.INIT_TAKEOFF = 0
            end
        end

        XLSX.openxlsx(file_xlsx, mode="rw") do xf
            sheet = xf["Aircrafts"]
            
            # Effacer le contenu existant de la feuille
            XLSX.rename!(sheet, "Aircrafts")
            
            # Écrire le DataFrame
            XLSX.writetable!(sheet, collect(eachcol(df_ac)), names(df_ac))
        end

        file_json = folder_path_json*inst_name*".json"
        println("Traitement de: $file_json")
        data = JSON.parsefile(file_json)
        data_dict = Dict(data)
        
        for ac in data_dict["aircrafts"]
            if data_dict["initial_flying_time"][ac] == 0
                data_dict["initial_takeoff"][ac] = 0
            end 
        end        
        # 4. Sauvegarder
        open(file_json, "w") do f
            JSON.print(f, data_dict;)
        end

    end  
 
end