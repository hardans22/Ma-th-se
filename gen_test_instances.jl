using CSV, DataFrames, Dates, JSON, XLSX, StatsBase


function detection(df_flights, aircrafts) 
    ac_problem = []
    for ac in aircrafts
        ac_df = filter(row -> row.TAIL_NUMBER == ac, df_flights)
        ac_df = sort(ac_df, :DEPARTURE_TIME)
        #println(ac_df)
        ac_list = [ac_df[1,"ORIGIN_AIRPORT"]]
        for i in 2:size(ac_df,1)
            if ac_df[i-1, "DESTINATION_AIRPORT"] == ac_df[i, "ORIGIN_AIRPORT"] &&  35 <= ac_df[i, "DEPARTURE_TIME"]-ac_df[i-1, "ARRIVAL_TIME"] && ac_df[i, "DEPARTURE_TIME"]-ac_df[i-1, "ARRIVAL_TIME"] <= 1440  
                push!(ac_list, ac_df[i-1, "DESTINATION_AIRPORT"])
            else
                push!(ac_problem,ac)
            end
        end
    end 
    return unique(ac_problem)
end

function process_xlsx_files(folder_path)
    # Obtenir tous les fichiers .xlsx
    xlsx_files = filter(f -> endswith(f, ".xlsx"), readdir(folder_path))
    
    for filename in xlsx_files
        filepath = joinpath(folder_path, filename)
        #println("Traitement de: $filename")
        
        # Traiter le fichier
        # wb = XLSX.readxlsx(filepath)
        # ... votre code de traitement ...
    end
    return xlsx_files
end

function process_xlsx_files_timer(folder_path)
    # Obtenir tous les fichiers .xlsx
    xlsx_files = filter(f -> endswith(f, ".xlsx"), readdir(folder_path))
    
    # Temps actuel
    current_time = time()
    
    # Filtrer les fichiers des 2 dernières minutes (120 secondes)
    recent_files = filter(xlsx_files) do filename
        filepath = joinpath(folder_path, filename)
        file_mtime = mtime(filepath)  # Temps de modification
        (current_time - file_mtime) <= 1200 # 120 secondes = 2 minutes
    end
    
    for filename in recent_files
        filepath = joinpath(folder_path, filename)
        println("Traitement de: $filename")
        # Traiter le fichier
        # wb = XLSX.readxlsx(filepath)
        # ... votre code de traitement ...
    end
    
    return recent_files
end

function computation(df_flights)
    tail_numbers = unique(df_flights.TAIL_NUMBER)
    aircraft_day_df = combine(groupby(df_flights, [:TAIL_NUMBER, :DAY]),
        :AIR_TIME => sum => :FLYING_TIME,
        :AIR_TIME => length => :TAKEOFF
    )

    sort!(aircraft_day_df, [:TAIL_NUMBER, :DAY])
    nbr_ac = length(tail_numbers)
    fh_ac_day = round(Int, sum(aircraft_day_df.FLYING_TIME)/nbr_ac/7)
    tk_ac_day = round(Int, sum(aircraft_day.TAKEOFF)/nbr_ac/7)
    fh_tk = round(sum(Int, aircraft_day.FLYING_TIME)/sum(aircraft_day.TAKEOFF))

    return (fh_ac_day = fh_ac_day, tk_ac_day = tk_ac_day, fh_tk = fh_tk)
end 

function computation_new(df_flights, nbr_ac)
    # Agrégation par jour pour obtenir les totaux quotidiens
    daily_stats = combine(groupby(df_flights, :DAY),
        :AIR_TIME => sum => :FLYING_TIME,
        :AIR_TIME => length => :TAKEOFF
    )
    
    # Calcul des moyennes par avion et par jour
    total_flying_time = sum(daily_stats.FLYING_TIME)
    total_takeoffs = sum(daily_stats.TAKEOFF)
    nbr_days = length(unique(df_flights.DAY))
    
    fh_ac_day = round(Int, total_flying_time / nbr_ac / nbr_days)
    tk_ac_day = round(Int, total_takeoffs / nbr_ac / nbr_days)
    fh_tk = round(Int, total_flying_time / total_takeoffs, digits=2)
    
    return (fh_ac_day = fh_ac_day, tk_ac_day = tk_ac_day, fh_tk = fh_tk)
end


#= month = "01"
file = "AS_2024-"*month*"_real_fl"
df_flights = DataFrame(XLSX.readtable("./Notebook/" * file * ".xlsx", "Sheet1"))

first_fold = "INSTANCES/instances_xlsx/"

tail_numbers = unique(df_flights.TAIL_NUMBER)
for nbr_day in [1]
    week = Int((nbr_day+6)/7)
    # ------------------- Part 1 -------------------
    list_ac = []
    for ac in tail_numbers
        ac_df = filter(row -> row.TAIL_NUMBER == ac, df_flights)
        ac_days = unique(ac_df.DAY)
        if issubset(collect(nbr_day:nbr_day+6), ac_days)
            push!(list_ac, ac)
        end
    end

    df_fl = filter(row -> (nbr_day <= row.DAY <= nbr_day+6) && (row.TAIL_NUMBER in list_ac), df_flights)
    ac_problem = detection(df_fl, list_ac)  # tu dois définir cette fonction ailleurs
    good_ac = setdiff(list_ac, ac_problem)
    # ------------------- Part 2 -------------------

    for len_ac in [1]
        aircrafts = good_ac[2:1+len_ac]
        #aircrafts = sample(good_ac, len_ac, replace=false)
        df_fl_ac = filter(row -> row.TAIL_NUMBER in aircrafts && row.DAY in nbr_day:nbr_day+6, df_flights)
        nbr_flights = size(df_fl_ac, 1)

        origin_airports = df_fl_ac.ORIGIN_AIRPORT
        sorted_o_apt = sort(collect(countmap(origin_airports)), by = x -> x[2], rev = true)
        cap_mat = [rand() < 0.2 ? 0 : rand(1:3) for _ in 1:length(sorted_o_apt), _ in 1:7]
        nbr_mtn_st = 7
        mtn_stations = [x[1] for x in sorted_o_apt[1:nbr_mtn_st]]

        # Générer capacité par jour
        m_st_df = DataFrame(MTN_STATIONS = mtn_stations)
        cap_bis = cap_mat[1:nbr_mtn_st,:]
        for i in 1:7
            m_st_df[!, "T_" * string(i)] = cap_bis[:, i]
        end

        fl_ind = computation(df_flights)
        fh_day = fl_ind.fh_ac_day
        fh_tk = fl_ind.fh_tk
        tk_day =fl_ind.tk_ac_day
        
        for version in 1:20
            for ac_critique in [1]
                init_apt, init_fl, init_tk, init_fd = [], [], [], []
                critical_indices = sample(1:length(aircrafts), ac_critique, replace=false)
                for (i, ac) in enumerate(aircrafts)
                    # Assigner rand(1800:4200) si l'avion est dans la sélection critique, sinon 0
                    temp_fl = (i in critical_indices) ? rand(2400:4200) : 0
                    
                    push!(init_apt, filter(row -> row.TAIL_NUMBER == ac, df_fl_ac)[1, "ORIGIN_AIRPORT"])
                    push!(init_fl, temp_fl)
                    push!(init_tk, temp_fl == 0 ? 1 : ceil(Int, temp_fl / 120)) # 2h (120 min) en moyenne par vol
                    push!(init_fd, temp_fl == 0 ? 1 : ceil(Int, temp_fl / 600)) # 10h (600 min) en moyenne par jour
                end
            
                Output_fold = first_fold*"A_MTN_"*string(ac_critique)*"/"
                # ------------------- Part 3 -------------------
                prm_df = DataFrame(TRT = 35, F = 6000, MT = 480, T = 50, D = 10, NBR_TP = 7, FH_TK = fh_tk, FH_DAY = fh_day, TK_DAY = tk_day)
                aircraft_df = DataFrame(
                    TAIL_NUMBER = aircrafts,
                    INIT_AIRPORT = init_apt,
                    INIT_FLYING_TIME = init_fl,
                    INIT_TAKEOFF = init_tk,
                    INIT_FLYING_DAY = init_fd
                )

                filename = Output_fold*string(nbr_flights)*"FL_"*string(len_ac)*"A_"*string(version)*".xlsx"
                XLSX.openxlsx(filename, mode = "w") do xf
                    # Supprimer la feuille par défaut "Sheet1"
                    XLSX.rename!(xf["Sheet1"], "Data")
                    data_sheet = xf["Data"]
        
                    #data_sheet = XLSX.addsheet!(xf, "Data")
                    parameters_sheet = XLSX.addsheet!(xf, "Parameters")
                    mtn_stations_sheet = XLSX.addsheet!(xf, "M_stations")
                    aircrafts_sheet = XLSX.addsheet!(xf, "Aircrafts")

                    XLSX.writetable!(data_sheet, Tables.columntable(df_fl_ac); write_columnnames = true)
                    XLSX.writetable!(parameters_sheet, Tables.columntable(prm_df); write_columnnames = true)
                    XLSX.writetable!(mtn_stations_sheet, Tables.columntable(m_st_df); write_columnnames = true)
                    XLSX.writetable!(aircrafts_sheet, Tables.columntable(aircraft_df); write_columnnames = true)
                    println("Fichier xlsx créé")
                end
            end
        end
    end
end
 =#
println()
for ac_critique in [15]
    folder_path ="INSTANCES/instances_literature_xlsx/A_MTN_"*string(ac_critique)*"/"  # ou le chemin vers votre dossier
    xlsx_files = process_xlsx_files_timer(folder_path)

    for file in xlsx_files
        df_flight = DataFrame(XLSX.readtable(folder_path*file, "Data"))
        df_param = DataFrame(XLSX.readtable(folder_path*file, "Parameters"))
        df_mstations = DataFrame(XLSX.readtable(folder_path*file, "M_stations"))
        df_aircrafts = DataFrame(XLSX.readtable(folder_path*file, "Aircrafts"))
        #df_inventory = DataFrame(XLSX.readtable(folder_path*file*".xlsx", "Inventory"))

        O_airport = unique(df_flight.ORIGIN_AIRPORT)
        D_airport = unique(df_flight.DESTINATION_AIRPORT)
        airports = unique(vcat(O_airport, D_airport))
        aircrafts = unique(df_aircrafts.TAIL_NUMBER)
        nbr_flights = nrow(df_flight)
        nbr_airports = length(airports)
        nbr_aircrafts = length(aircrafts)

        mtn_stations = df_mstations.MTN_STATIONS
        initial_airport = Dict(row.TAIL_NUMBER => row.INIT_AIRPORT for row in eachrow(df_aircrafts))
        nbr_mstations = length(mtn_stations)
        initial_flying_time = Dict(row.TAIL_NUMBER => row.INIT_FLYING_TIME for row in eachrow(df_aircrafts))
        initial_takeoff = Dict(row.TAIL_NUMBER => row.INIT_TAKEOFF for row in eachrow(df_aircrafts))
        initial_flying_day = Dict(row.TAIL_NUMBER => row.INIT_FLYING_DAY for row in eachrow(df_aircrafts))
        #ms_capacity = Dict(row.MTN_STATIONS => [] for row in eachrow(df_mstations))
        ms_capacity = Dict(row.MTN_STATIONS => collect(row[2:end]) for row in eachrow(df_mstations))
        turn_around_time = df_param.TRT[1]
        flying_time_max = df_param.F[1]
        takeoff_max = df_param.T[1]
        flying_day_max = df_param.D[1]
        mtn_time = df_param.MT[1]
        nbr_TP = df_param.NBR_TP[1]
        fh_tk = df_param.FH_TK[1]
        fh_day = df_param.FH_DAY[1]
        tk_day = df_param.TK_DAY[1]

        #= 
        exp_part = df_inventory.EXP_PART
        nbr_exp_part = length(exp_part)
        init_level_ep = Dict(row.EXP_PART => row.INIT_LEVEL for row in eachrow(df_inventory))
        rate_ep = Dict(row.EXP_PART => row.RATE for row in eachrow(df_inventory))
        =#
        df_copy = deepcopy(df_flight)
        #df_copy = time_to_minutes(df_copy)
        temp1 = maximum(df_flight.DAY)*1440
        temp2 = maximum(df_copy.ARRIVAL_TIME)
        end_horizon_time = maximum([temp1, temp2])

        flight_legs = []
        for row in eachrow(df_copy)
            push!(flight_legs, (String(row.ORIGIN_AIRPORT), String(row.DESTINATION_AIRPORT), Int(row.DEPARTURE_TIME), Int(row.ARRIVAL_TIME), Int(row.AIR_TIME)))
        end

        instance_data = Dict(
            "number_of_flight_legs" => nbr_flights, "number_of_airports" => nbr_airports,
            "number_of_maintenance_stations" => nbr_mstations, "number_of_aircrafts" => nbr_aircrafts,
            "aircrafts" => aircrafts, "maximum_flying_time" => flying_time_max,
            "maximum_takeoff" => takeoff_max, "maximum_flying_day" => flying_day_max, "airports" => airports,
            "maintenance_stations" => mtn_stations, "initial_flying_time" => initial_flying_time,
            "initial_takeoff" => initial_takeoff, "initial_flying_day" => initial_flying_day,
            "mtn_station_capacity" => ms_capacity, "initial_airport_aircraft" => initial_airport, 
            "flight_legs" => flight_legs, "turn_around_time" => turn_around_time, 
            "maintenance_time" => mtn_time, "end_horizon_time" => end_horizon_time,
            "nbr_TP" => nbr_TP, "fh_tk" => fh_tk, "fh_day" => fh_day, "tk_day" => tk_day
            #= "exp_part" => exp_part,
            "number_of_exp_part" => nbr_exp_part,
            "init_level_ep" => init_level_ep,
            "rate_ep" => rate_ep =#
        )

        # S'assurer que le dossier existe avant d'écrire le fichier JSON
        json_filepath = "INSTANCES/instances_literature_json/A_MTN_"*string(ac_critique)*"/" * string(splitext(file)[1]) * ".json"
        dirpath = dirname(json_filepath)
        if !isdir(dirpath)
            mkpath(dirpath)
        end

        # Sauvegarde au format JSON dans un fichier
        open(json_filepath, "w") do f
            JSON.print(f, instance_data;)
            println("Fichier json créé")
        end
    end
end