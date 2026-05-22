using JuMP, Gurobi, JSON, MathOptInterface, DataFrames, XLSX, CSV, Dates

include("../../build_graph.jl")
include("__benders__.jl")
include("__cuts_func__.jl")

env = Gurobi.Env()

instance = ARGS[1]
sp_method = ARGS[2]

nbr_thread = 14
silent = false
time_limit = 1800
cut_type = "desagg"

# =============== Sur Nibi ===============
#instance = ARGS[1]
#sp_method = ARGS[2]


# =============== En local sur mon ordi ===============
inst = ARGS[1]
sp_method = ARGS[2]
instance = "../../INSTANCES/instances_new_json/"*inst

# ============== Gestion des fichiers et la sortie =================
parts = split(instance, "/")
f_fold = parts[end-1]
Outputs_fold = "RESULTS_BENDERS/" * f_fold * "/"
inst_name = splitext(basename(instance))[1]


path_file = Outputs_fold*"result_"*inst_name*"_benders_"*sp_method*".txt"
Output_file = open(path_file, "w") 
write_both(Output_file, "-----------------------INSTANCE $inst_name--------------------------")
close(Output_file)

instance_data = build_graph(instance)
solution = benders_decomp(env, instance_data, sp_method, cut_type, path_file, nbr_thread, silent, time_limit)
other_info = solution.other_info
Output_file = open(path_file, "a") 
#print_solution(solution, instance_data, Output_file)

Obj = solution.obj
MP_time = other_info["mp_time"]
SP_time = other_info["sp_time"]
Time = other_info["time"]
nbr_feas_cuts = other_info["nbr_feas_cuts"]
nbr_opti_cuts = other_info["nbr_opti_cuts"]
nbr_fix_cuts = other_info["nbr_fix_cuts"]
nbr_cuts = other_info["nbr_cuts"]
nbr_iter = other_info["nbr_iter"]
if solution.other_info["status"] != "NO SOLUTION FOUND"
    nbr_mtn = other_info["nbr_mtn"]
else
    nbr_mtn = 0
end

write_both(Output_file, "\n=== Résultat final ===")
write_both(Output_file, "Objectif: $(solution.obj)")
write_both(Output_file, "Temps total du master: $(round(other_info["mp_time"], digits=4))s")
write_both(Output_file, "Temps total du sp: $(round(other_info["sp_time"], digits=4))s") 
write_both(Output_file, "Temps total: $(other_info["time"])s")
write_both(Output_file, "Nombre de coupes: $(other_info["nbr_cuts"])")
write_both(Output_file, "Nombre de coupes faisabilité: $(other_info["nbr_feas_cuts"])")
write_both(Output_file, "Nombre de coupes d'optimalité: $(other_info["nbr_opti_cuts"])")
write_both(Output_file, "Nombre de coupes de fixing: $(other_info["nbr_fix_cuts"])")
write_both(Output_file, "Nombre d'itérations: $(other_info["nbr_iter"])")

dataframe = DataFrames.DataFrame(Instances = [inst_name], UB = [Obj], MP_Time = [MP_time], SP_Time = [SP_time],
                                Time = [Time], Nbr_mtn = nbr_mtn, Nbr_iter = [nbr_iter], Nbr_feas_cuts = [nbr_feas_cuts],
                                Nbr_opti_cuts = [nbr_opti_cuts],  Nbr_fix_cuts = nbr_fix_cuts, Nbr_cuts = [nbr_cuts])

csv_summary_path = Outputs_fold*"summary_"*inst_name*"_benders_"*sp_method*".csv"
# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

println("CSV Summary saved to: ", csv_summary_path)
