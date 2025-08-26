using JuMP
using Gurobi
using JSON
using MathOptInterface
using DataFrames
using XLSX

include("model_amrp.jl")
inst_name = "63FL_5A_3D_3"
file_path = "instances_json/"*inst_name
nbr_thread = 10
silent = false
time_limit = 150

println("INSTANCE : $file_path")
instance = JSON.parsefile(file_path*".json")
compt = 0
const MOI = MathOptInterface
while compt < 1
    global compt 
    new_instance = deepcopy(instance)
    # Modifier les valeurs dans initial_flying_time
    for key in keys(new_instance["initial_flying_time"])
        new_instance["initial_flying_time"][key] = rand(1500:3000)
        new_instance["initial_takeoff"][key] = Int(ceil(new_instance["initial_flying_time"][key]/180))
        new_instance["initial_flying_day"][key] = Int(ceil(new_instance["initial_flying_time"][key]/600)) 
    end

    inst_dict = build_graph(new_instance)
    solution = model_amrp(inst_dict, nbr_thread, silent, time_limit)
 
    if solution["status"] == MOI.OPTIMAL || solution["status"] == MOI.TIME_LIMIT || solution["status"] == MOI.FEASIBLE_POINT
        println("-------------------UNE INSTANCE TROUVÉE ---------------------")
        compt+=1
        # Écrire la nouvelle instance dans un nouveau fichier
        open(file_path * "_" * string(compt) * ".json", "w") do f
            JSON.print(f, new_instance)
        end
    end 

end