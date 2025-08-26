using CSV
using DataFrames
using Dates
using JSON
using XLSX

function time_to_minutes(df)
    for row in eachrow(df_copy)
        #println(row.DEPARTURE_TIME)
        dt = hour(Time(row.DEPARTURE_TIME, "H:M"))*60 + 1440*(row.DAY-1) + minute(Time(row.DEPARTURE_TIME, "H:M"))
        at = hour(Time(row.ARRIVAL_TIME, "H:M"))*60 + 1440*(row.DAY-1) + minute(Time(row.ARRIVAL_TIME, "H:M"))
         
        #= dt = hour(row.DEPARTURE_TIME)*60 + 1440*(row.DAY-1) + minute(row.DEPARTURE_TIME)
        at = hour(row.ARRIVAL_TIME)*60 + 1440*(row.DAY-1) + minute(row.ARRIVAL_TIME)
         =#
        if at < dt
            at += 1440
        end 
        row.DEPARTURE_TIME = dt
        row.ARRIVAL_TIME = at
    end
    return df 
end
  
#file = "47FL_5A_3D" 
for f in ["171FL_16A_3D_"]
    for version in 1:3
        file = f*string(version)
        df_flights = DataFrame(XLSX.readtable("instances_xlsx/"*file*".xlsx", "Data"))
        df_param = DataFrame(XLSX.readtable("instances_xlsx/"*file*".xlsx", "Parameters"))
        df_mstations = DataFrame(XLSX.readtable("instances_xlsx/"*file*".xlsx", "M_stations"))
        df_aircrafts = DataFrame(XLSX.readtable("instances_xlsx/"*file*".xlsx", "Aircrafts"))
        #df_inventory = DataFrame(XLSX.readtable("instances_xlsx/"*file*".xlsx", "Inventory"))

        O_airport = unique(df_flights.ORIGIN_AIRPORT)
        D_airport = unique(df_flights.DESTINATION_AIRPORT)
        airports = unique(vcat(O_airport, D_airport))
        aircrafts = unique(df_flights.TAIL_NUMBER)
        nbr_flights = nrow(df_flights)
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
        df_copy = deepcopy(df_flights)
        #df_copy = time_to_minutes(df_copy)
        temp1 = maximum(df_flights.DAY)*1440
        temp2 = maximum(df_copy.ARRIVAL_TIME)
        end_horizon_time = maximum([temp1, temp2])

        flight_legs = []
        for row in eachrow(df_copy)
            push!(flight_legs, (String(row.ORIGIN_AIRPORT), String(row.DESTINATION_AIRPORT), Int(row.DEPARTURE_TIME), Int(row.ARRIVAL_TIME), Int(row.AIR_TIME)))
        end

        instance_data = Dict(
            "number_of_flight_legs" => nbr_flights,
            "number_of_airports" => nbr_airports,
            "number_of_maintenance_stations" => nbr_mstations,
            "number_of_aircrafts" => nbr_aircrafts,
            "aircrafts" => aircrafts,
            "maximum_flying_time" => flying_time_max,
            "maximum_takeoff" => takeoff_max,
            "maximum_flying_day" => flying_day_max,
            "airports" => airports,
            "maintenance_stations" => mtn_stations,
            "initial_flying_time" => initial_flying_time,
            "initial_takeoff" => initial_takeoff,
            "initial_flying_day" => initial_flying_day,
            "mtn_station_capacity" => ms_capacity,
            "initial_airport_aircraft" => initial_airport, 
            "flight_legs" => flight_legs,
            "turn_around_time" => turn_around_time,
            "maintenance_time" => mtn_time,
            "end_horizon_time" => end_horizon_time,
            "nbr_TP" => nbr_TP
            #= "exp_part" => exp_part,
            "number_of_exp_part" => nbr_exp_part,
            "init_level_ep" => init_level_ep,
            "rate_ep" => rate_ep =#
        )

        # Sauvegarde au format JSON dans un fichier
        open("instances_json/"*file*".json", "w") do f
            JSON.print(f, instance_data;)
        end
    end
end