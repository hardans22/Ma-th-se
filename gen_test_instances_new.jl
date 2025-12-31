using CSV, DataFrames, Dates, JSON, XLSX, StatsBase

# ============= Fonctions utilitaires =============

function detect_problematic_aircraft(df_flights, aircrafts)
    """Détecte les avions avec des problèmes de continuité de vol"""
    problematic = Set{String}()
    
    for ac in aircrafts
        ac_df = sort(filter(row -> row.TAIL_NUMBER == ac, df_flights), :DEPARTURE_TIME)
        
        for i in 2:nrow(ac_df)
            turnaround = ac_df[i, :DEPARTURE_TIME] - ac_df[i-1, :ARRIVAL_TIME]
            different_airport = ac_df[i-1, :DESTINATION_AIRPORT] != ac_df[i, :ORIGIN_AIRPORT]
            invalid_turnaround = !(35 <= turnaround <= 1440)
            
            if different_airport || invalid_turnaround
                push!(problematic, ac)
                break
            end
        end
    end
    
    return collect(problematic)
end

function get_recent_xlsx_files(folder_path, max_age_seconds=120)
    """Retourne les fichiers xlsx modifiés récemment"""
    current_time = time()
    xlsx_files = filter(f -> endswith(f, ".xlsx"), readdir(folder_path))
    
    return filter(xlsx_files) do filename
        filepath = joinpath(folder_path, filename)
        (current_time - mtime(filepath)) <= max_age_seconds
    end
end

function compute_flight_indicators(df_flights)
    """Calcule les indicateurs de vol moyens"""
    aircraft_day_df = combine(
        groupby(df_flights, [:TAIL_NUMBER, :DAY]),
        :AIR_TIME => sum => :FLYING_TIME,
        :AIR_TIME => length => :TAKEOFF
    )
    
    nbr_ac = length(unique(df_flights.TAIL_NUMBER))
    nbr_days = 7
    
    return (
        fh_ac_day = round(Int, sum(aircraft_day_df.FLYING_TIME) / nbr_ac / nbr_days),
        tk_ac_day = round(Int, sum(aircraft_day_df.TAKEOFF) / nbr_ac / nbr_days),
        fh_tk = round(Int, sum(aircraft_day_df.FLYING_TIME) / sum(aircraft_day_df.TAKEOFF))
    )
end

function filter_aircraft_by_week(df_flights, tail_numbers, start_day)
    """Retourne les avions disponibles toute la semaine"""
    week_days = start_day:start_day+6
    
    return filter(tail_numbers) do ac
        ac_days = unique(filter(row -> row.TAIL_NUMBER == ac, df_flights).DAY)
        issubset(week_days, ac_days)
    end
end

function generate_maintenance_stations(df_flights, nbr_stations=7)
    """Génère les stations de maintenance avec capacités"""
    origin_counts = sort(collect(countmap(df_flights.ORIGIN_AIRPORT)), by=x->x[2], rev=true)
    stations = [x[1] for x in origin_counts[1:nbr_stations]]
    
    # Capacité aléatoire par jour (0 ou 1-3)
    capacity_matrix = [rand() < 0.2 ? 0 : rand(1:3) for _ in 1:nbr_stations, _ in 1:7]
    
    df_stations = DataFrame(MTN_STATIONS = stations)
    for day in 1:7
        df_stations[!, "T_$day"] = capacity_matrix[:, day]
    end
    
    return df_stations
end

function generate_aircraft_initial_states(aircrafts, df_flights, nbr_critical=1)
    """Génère les états initiaux des avions"""
    critical_indices = Set(sample(1:length(aircrafts), nbr_critical, replace=false))
    
    init_data = map(enumerate(aircrafts)) do (i, ac)
        is_critical = i in critical_indices
        flying_time = is_critical ? rand(2400:4200) : 0
        
        (
            TAIL_NUMBER = ac,
            INIT_AIRPORT = filter(row -> row.TAIL_NUMBER == ac, df_flights)[1, :ORIGIN_AIRPORT],
            INIT_FLYING_TIME = flying_time,
            INIT_TAKEOFF = flying_time == 0 ? 1 : ceil(Int, flying_time / 120),
            INIT_FLYING_DAY = flying_time == 0 ? 1 : ceil(Int, flying_time / 600)
        )
    end
    
    return DataFrame(init_data)
end

function create_instance_xlsx(filename, df_flights, df_params, df_stations, df_aircraft)
    """Crée un fichier xlsx avec toutes les feuilles nécessaires"""
    XLSX.openxlsx(filename, mode="w") do xf
        XLSX.rename!(xf["Sheet1"], "Data")
        
        sheets = Dict(
            "Data" => df_flights,
            "Parameters" => df_params,
            "M_stations" => df_stations,
            "Aircrafts" => df_aircraft
        )
        
        for (name, df) in sheets
            sheet = name == "Data" ? xf["Data"] : XLSX.addsheet!(xf, name)
            XLSX.writetable!(sheet, Tables.columntable(df), write_columnnames=true)
        end
    end
end

function xlsx_to_json(xlsx_path, json_path)
    """Convertit un fichier xlsx en format json"""
    # Charger les données
    df_flight = DataFrame(XLSX.readtable(xlsx_path, "Data"))
    df_param = DataFrame(XLSX.readtable(xlsx_path, "Parameters"))
    df_mstations = DataFrame(XLSX.readtable(xlsx_path, "M_stations"))
    df_aircrafts = DataFrame(XLSX.readtable(xlsx_path, "Aircrafts"))
    
    # Extraire les informations
    airports = unique(vcat(df_flight.ORIGIN_AIRPORT, df_flight.DESTINATION_AIRPORT))
    aircrafts = df_aircrafts.TAIL_NUMBER
    mtn_stations = df_mstations.MTN_STATIONS
    
    # Créer les dictionnaires
    initial_airport = Dict(row.TAIL_NUMBER => row.INIT_AIRPORT for row in eachrow(df_aircrafts))
    initial_flying_time = Dict(row.TAIL_NUMBER => row.INIT_FLYING_TIME for row in eachrow(df_aircrafts))
    initial_takeoff = Dict(row.TAIL_NUMBER => row.INIT_TAKEOFF for row in eachrow(df_aircrafts))
    initial_flying_day = Dict(row.TAIL_NUMBER => row.INIT_FLYING_DAY for row in eachrow(df_aircrafts))
    ms_capacity = Dict(row.MTN_STATIONS => collect(row[2:end]) for row in eachrow(df_mstations))
    
    # Créer les flight legs
    flight_legs = [(
        String(row.ORIGIN_AIRPORT), 
        String(row.DESTINATION_AIRPORT), 
        Int(row.DEPARTURE_TIME), 
        Int(row.ARRIVAL_TIME), 
        Int(row.AIR_TIME)
    ) for row in eachrow(df_flight)]
    
    # Calculer l'horizon temporel
    end_horizon = max(maximum(df_flight.DAY) * 1440, maximum(df_flight.ARRIVAL_TIME))
    
    # Construire le dictionnaire d'instance
    instance_data = Dict(
        "number_of_flight_legs" => nrow(df_flight),
        "number_of_airports" => length(airports),
        "number_of_maintenance_stations" => length(mtn_stations),
        "number_of_aircrafts" => length(aircrafts),
        "aircrafts" => aircrafts,
        "airports" => airports,
        "maintenance_stations" => mtn_stations,
        "maximum_flying_time" => df_param.F[1],
        "maximum_takeoff" => df_param.T[1],
        "maximum_flying_day" => df_param.D[1],
        "initial_flying_time" => initial_flying_time,
        "initial_takeoff" => initial_takeoff,
        "initial_flying_day" => initial_flying_day,
        "mtn_station_capacity" => ms_capacity,
        "initial_airport_aircraft" => initial_airport,
        "flight_legs" => flight_legs,
        "turn_around_time" => df_param.TRT[1],
        "maintenance_time" => df_param.MT[1],
        "end_horizon_time" => end_horizon,
        "nbr_TP" => df_param.NBR_TP[1],
        "fh_tk" => df_param.FH_TK[1],
        "fh_day" => df_param.FH_DAY[1],
        "tk_day" => df_param.TK_DAY[1]
    )
    
    # Sauvegarder en JSON
    mkpath(dirname(json_path))
    open(json_path, "w") do f
        JSON.print(f, instance_data)
    end
end

# ============= Script principal =============

function main()
    month = "01"
    file = "AS_2024-$(month)_real_fl"
    df_flights = DataFrame(XLSX.readtable("./Notebook/$(file).xlsx", "Sheet1"))
    
    base_folder = "INSTANCES/instances_xlsx/"
    tail_numbers = unique(df_flights.TAIL_NUMBER)
    
    # Paramètres de génération
    nbr_days_list = [1]
    nbr_aircraft_list = [1]
    nbr_versions = 20
    nbr_critical_aircraft_list = [1]
    
    # Indicateurs de vol
    flight_indicators = compute_flight_indicators(df_flights)
    
    for start_day in nbr_days_list
        # Filtrer les avions disponibles toute la semaine
        available_aircraft = filter_aircraft_by_week(df_flights, tail_numbers, start_day)
        df_week = filter(row -> (start_day <= row.DAY <= start_day+6) && (row.TAIL_NUMBER in available_aircraft), df_flights)
        
        # Détecter les avions problématiques
        problematic = detect_problematic_aircraft(df_week, available_aircraft)
        good_aircraft = setdiff(available_aircraft, problematic)
        
        for nbr_aircraft in nbr_aircraft_list
            # Sélectionner les avions
            selected_aircraft = good_aircraft[1:min(nbr_aircraft, length(good_aircraft))]
            df_selected = filter(row -> row.TAIL_NUMBER in selected_aircraft, df_week)
            
            # Générer les stations de maintenance
            df_stations = generate_maintenance_stations(df_selected)
            
            for nbr_critical in nbr_critical_aircraft_list
                output_folder = "$(base_folder)A_MTN_$(nbr_critical)/"
                mkpath(output_folder)
                
                for version in 1:nbr_versions
                    # Générer les états initiaux des avions
                    df_aircraft = generate_aircraft_initial_states(selected_aircraft, df_selected, nbr_critical)
                    
                    # Créer le DataFrame de paramètres
                    df_params = DataFrame(
                        TRT = 35, F = 6000, MT = 480, T = 50, D = 10, 
                        NBR_TP = 7, FH_TK = flight_indicators.fh_tk,
                        FH_DAY = flight_indicators.fh_ac_day, TK_DAY = flight_indicators.tk_ac_day
                    )
                    
                    # Sauvegarder en xlsx
                    filename = "$(output_folder)$(nrow(df_selected))FL_$(nbr_aircraft)A_$(version).xlsx"
                    create_instance_xlsx(filename, df_selected, df_params, df_stations, df_aircraft)
                    println("Fichier xlsx créé: $filename")
                end
            end
        end
    end
    
    # Conversion xlsx -> json
    println("\n=== Conversion xlsx vers json ===")
    for nbr_critical in nbr_critical_aircraft_list
        xlsx_folder = "$(base_folder)A_MTN_$(nbr_critical)/"
        json_folder = "INSTANCES/instances_json/A_MTN_$(nbr_critical)/"
        
        recent_files = get_recent_xlsx_files(xlsx_folder)
        
        for file in recent_files
            xlsx_path = joinpath(xlsx_folder, file)
            json_path = joinpath(json_folder, replace(file, ".xlsx" => ".json"))
            
            xlsx_to_json(xlsx_path, json_path)
            println("Fichier json créé: $json_path")
        end
    end
end

# Exécuter le script
main()