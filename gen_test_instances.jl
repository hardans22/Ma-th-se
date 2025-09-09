using CSV
using DataFrames
using Dates
using JSON
using XLSX
using StatsBase


function detection(df_flights, aircrafts) 
    ac_problem = []
    for ac in aircrafts
        ac_df = filter(row -> row.TAIL_NUMBER == ac, df_flights)
        #println(ac_df)
        ac_list = [ac_df[1,"ORIGIN_AIRPORT"]]
        for i in 2:size(ac_df,1)
            if ac_df[i-1, "DESTINATION_AIRPORT"] == ac_df[i, "ORIGIN_AIRPORT"] &&  30 <= ac_df[i, "DEPARTURE_TIME"]-ac_df[i-1, "ARRIVAL_TIME"] && ac_df[i, "DEPARTURE_TIME"]-ac_df[i-1, "ARRIVAL_TIME"] <= 1440  
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


month = "01"
file = "AS_2015-"*month*"_real_fl"
df_flights = DataFrame(XLSX.readtable("./Notebook/" * file * ".xlsx", "Sheet1"))

first_fold = "instances_xlsx/M"*month

tail_numbers = unique(df_flights.TAIL_NUMBER)
for nbr_day in [1,8,15,22]
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

    df_fl = filter(row -> row.TAIL_NUMBER in list_ac, df_flights)
    ac_problem = detection(df_fl, list_ac)  # tu dois définir cette fonction ailleurs
    good_ac = setdiff(list_ac, ac_problem)
    len_ac = length(good_ac)

    # ------------------- Part 2 -------------------
    
    for len in [5,10,15]
        aircrafts = good_ac[1:len]
        #aircrafts = sample(good_ac, len, replace=false)
        df_fl_ac = filter(row -> row.TAIL_NUMBER in aircrafts && row.DAY in nbr_day:nbr_day+6, df_flights)
        nbr_flights = size(df_fl_ac, 1)
        nbr_aircrafts = length(aircrafts)

        origin_airports = df_fl_ac.ORIGIN_AIRPORT
        sorted_o_apt = sort(collect(countmap(origin_airports)), by = x -> x[2], rev = true)
        cap_mat = [rand() < 0.2 ? 0 : rand(1:3) for _ in 1:length(sorted_o_apt), _ in nbr_day:nbr_day+6]

        v_range = len <= 5 ? (1:2) : (1:3)
        for version in v_range
            init_apt, init_fl, init_tk, init_fd = [], [], [], []
            nbr_ac_mtn = [3, 5, 8][version]
            for (i, ac) in enumerate(aircrafts)
                temp_fl = (i <= length(aircrafts) - nbr_ac_mtn) ? 0 : rand(1800:4200)
                
                push!(init_apt, filter(row -> row.TAIL_NUMBER == ac, df_fl_ac)[1, "ORIGIN_AIRPORT"])
                push!(init_fl, temp_fl)
                push!(init_tk, temp_fl == 0 ? 1 : ceil(Int, temp_fl / 150))     # 2.5h (150 min) en moyenne par vol
                push!(init_fd, temp_fl == 0 ? 1 : ceil(Int, temp_fl / 600))     # 10h (600 min) en moyenne par jour
            end 
            for nbr_mtn_st in [5, 10]
                Output_fold = first_fold*"_W"*string(week)*"/MS_"*string(nbr_mtn_st)*"/"
                mtn_stations = [x[1] for x in sorted_o_apt[1:nbr_mtn_st]]

                # Générer capacité par jour
                m_st_df = DataFrame(MTN_STATIONS = mtn_stations)
                cap_bis = cap_mat[1:nbr_mtn_st,:]
                for i in 1:7
                    m_st_df[!, "T_" * string(i)] = cap_bis[:, i]
                end

                # ------------------- Part 3 -------------------
                prm_df = DataFrame(TRT = 30, F = 6000, MT = 480, T = 40, D = 10, NBR_TP = nbr_day+6)
                aircraft_df = DataFrame(
                    TAIL_NUMBER = aircrafts,
                    INIT_AIRPORT = init_apt,
                    INIT_FLYING_TIME = init_fl,
                    INIT_TAKEOFF = init_tk,
                    INIT_FLYING_DAY = init_fd
                )

                filename = Output_fold*string(nbr_flights)*"FL_"*string(nbr_aircrafts)*"A_"*string(version)*".xlsx"
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


println()
for nbr_day in [7,14,21,28]
    week = Int(nbr_day/7)
    for mtn_st in [5, 10]
        folder_path ="instances_xlsx/M"*month*"_W"*string(week)*"/MS_"*string(mtn_st)*"/"  # ou le chemin vers votre dossier
        xlsx_files = process_xlsx_files(folder_path)

        for file in xlsx_files
            df_flight = DataFrame(XLSX.readtable(folder_path*file, "Data"))
            df_param = DataFrame(XLSX.readtable(folder_path*file, "Parameters"))
            df_mstations = DataFrame(XLSX.readtable(folder_path*file, "M_stations"))
            df_aircrafts = DataFrame(XLSX.readtable(folder_path*file, "Aircrafts"))
            #df_inventory = DataFrame(XLSX.readtable(folder_path*file*".xlsx", "Inventory"))

            O_airport = unique(df_flight.ORIGIN_AIRPORT)
            D_airport = unique(df_flight.DESTINATION_AIRPORT)
            airports = unique(vcat(O_airport, D_airport))
            aircrafts = unique(df_flight.TAIL_NUMBER)
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
                "nbr_TP" => nbr_TP
                #= "exp_part" => exp_part,
                "number_of_exp_part" => nbr_exp_part,
                "init_level_ep" => init_level_ep,
                "rate_ep" => rate_ep =#
            )

            # Sauvegarde au format JSON dans un fichier
            open("instances_json/M"*month*"_W"*string(week)*"/MS_"*string(mtn_st)*"/"*splitext(file)[1]*".json", "w") do f
                JSON.print(f, instance_data;)
                println("Fichier json créé")
            end
        end
    end
end