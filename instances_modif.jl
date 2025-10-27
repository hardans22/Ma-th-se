using CSV, DataFrames, Dates, JSON, XLSX

function process_files(folder_path, file_type)
    # Obtenir tous les fichiers .xlsx
    files = filter(f -> endswith(f, file_type), readdir(folder_path))
    return files
end

function computation(df_flights)
    tail_numbers = unique(df_flights.TAIL_NUMBER)
    # air_time_by_tail_day est la somme des airtime de chaque tail pour chaque jour
    air_time_by_tail_day = Dict{Tuple{String, Int}, Float64}()

    # takeoff_by_tail_day est le nombre d'occurence de chaque tail number pour chaque jour
    takeoff_by_tail_day = Dict{Tuple{String, Int}, Int}()
    for row in eachrow(df_flights)
        tail = row.TAIL_NUMBER
        day = row.DAY
        airtime = row.AIR_TIME

        key = (tail, day)

        if haskey(air_time_by_tail_day, key)
            air_time_by_tail_day[key] += airtime
            takeoff_by_tail_day[key] += 1
        else
            air_time_by_tail_day[key] = airtime
            takeoff_by_tail_day[key] = 1
        end
    end

    result_df = DataFrame(TAIL_NUMBER = String[],DAY = Int[],TAKEOFF = Int[], FLYING_TIME = Float64[])

    for ((tail, day), takeoff) in sort(collect(takeoff_by_tail_day); by=x->(x[1][1], x[1][2]))
        flying_time = air_time_by_tail_day[(tail, day)]
        push!(result_df, (tail, day, takeoff, flying_time))
    end
    nbr_ac = length(tail_numbers)
    fh_ac_day = round(Int, sum(result_df.FLYING_TIME)/nbr_ac/7)
    tk_ac_day = round(Int, sum(result_df.TAKEOFF)/nbr_ac/7)
    fh_tk = round(sum(Int, result_df.FLYING_TIME)/sum(result_df.TAKEOFF))

    return Dict("fh_ac_day" => fh_ac_day, "tk_ac_day" => tk_ac_day, "fh_tk" => fh_tk)
end 

function computation2(df_flights, nbr_ac)
    # air_time_by_tail_day est la somme des airtime de chaque tail pour chaque jour
    air_time_by_tail_day = Dict{Tuple{String, Int}, Float64}()

    # Agrégation par jour uniquement
    takeoff_by_day = Dict{Int, Int}()
    air_time_by_day = Dict{Int, Float64}()

    for row in eachrow(df_flights)
        day = row.DAY
        airtime = row.AIR_TIME
        
        if haskey(air_time_by_day, day)
            air_time_by_day[day] += airtime
            takeoff_by_day[day] += 1
        else
            air_time_by_day[day] = airtime
            takeoff_by_day[day] = 1
        end
    end

    # Créer le DataFrame résultat
    result_df = DataFrame(DAY = Int[], TAKEOFF = Int[], FLYING_TIME = Float64[])

    for day in sort(collect(keys(takeoff_by_day)))
        takeoff = takeoff_by_day[day]
        flying_time = air_time_by_day[day]
        push!(result_df, (day, takeoff, flying_time))
    end

    fh_ac_day = round(Int, sum(result_df.FLYING_TIME)/nbr_ac/7)
    tk_ac_day = round(Int, sum(result_df.TAKEOFF)/nbr_ac/7)
    fh_tk = round(Int,sum(result_df.FLYING_TIME)/sum(result_df.TAKEOFF))

    return Dict("fh_ac_day" => fh_ac_day, "tk_ac_day" => tk_ac_day, "fh_tk" => fh_tk)
end 


for ac_critique in [2]  
    folder_path_xlsx ="INSTANCES/instances_xlsx/A_MTN_"*string(ac_critique)*"/"  # ou le chemin vers votre dossier
    folder_path_json ="INSTANCES/instances_json/A_MTN_"*string(ac_critique)*"/"  # ou le chemin vers votre dossier
    xlsx_files = process_files(folder_path_xlsx, "xlsx")
    for file in xlsx_files
        inst_name = split(file, ".")[1]
        file_xlsx = folder_path_xlsx*file
        println("Traitement de: $file_xlsx")
        df_flights = DataFrame(XLSX.readtable(file_xlsx, "Data"))
        df_param = DataFrame(XLSX.readtable(file_xlsx, "Parameters"))
        df_aircrafts = DataFrame(XLSX.readtable(file_xlsx, "Aircrafts"))
        
        nbr_ac = length(df_aircrafts.TAIL_NUMBER)
        df_param.AC_CRITIQUE .= ac_critique
        dict_gamma = computation(df_flights)
        df_param.FH_DAY .= dict_gamma["fh_ac_day"]
        df_param.FH_TK .= dict_gamma["fh_tk"]
        df_param.TK_DAY .= dict_gamma["tk_ac_day"]

        # Lire toutes les feuilles existantes
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

        file_json = folder_path_json*inst_name*".json"
        println("Traitement de: $file_json")
        data = JSON.parsefile(file_json)
        data_dict = Dict(data)
        data_dict["ac_critique"] = ac_critique
        data_dict["fh_day"] = dict_gamma["fh_ac_day"]
        data_dict["tk_day"] = dict_gamma["tk_ac_day"]
        data_dict["fh_tk"] = dict_gamma["fh_tk"]

        # 4. Sauvegarder
        open(file_json, "w") do f
            JSON.print(f, data_dict;)
        end

    end  
 
end