using JSON

nbr_instance = 10


for inst_name in ["146FL_15A_3D"]
    file = "instances_json/"*inst_name
    for i in 1:nbr_instance
        new_instance = JSON.parsefile(file*"_"*string(i)*".json")
        instance = JSON.parsefile(file*".json")
        # Modifier les valeurs dans initial_flying_time
        #= new_instance["nbr_TP"] = instance["nbr_TP"]
        new_instance["exp_part"] = instance["exp_part"]
        new_instance["number_of_exp_part"] = instance["number_of_exp_part"]
        new_instance["init_level_ep"] = instance["init_level_ep"]
        new_instance["rate_ep"] = instance["rate_ep"]
         =#
        # Écrire la nouvelle instance dans un nouveau fichier
        open(file * "_" * string(i) * ".json", "w") do f
            JSON.print(f, new_instance)
        end 

    end
end


#= for inst_name in ["47FL_5A_3D"]
    file = "instances_json/"*inst_name
    for i in 1:nbr_instance
        instance = JSON.parsefile(file*".json")
        new_instance = deepcopy(instance)
        # Modifier les valeurs dans initial_flying_time
        for key in keys(new_instance["initial_flying_time"])
            new_instance["initial_flying_time"][key] = rand(500:3000)
            new_instance["initial_takeoff"][key] = Int(ceil(new_instance["initial_flying_time"][key]/180))
            new_instance["initial_flying_day"][key] = Int(ceil(new_instance["initial_flying_time"][key]/600)) 
        end
    
        # Écrire la nouvelle instance dans un nouveau fichier
        open(file * "_" * string(i) * ".json", "w") do f
            JSON.print(f, new_instance)
        end 

    end
end =#