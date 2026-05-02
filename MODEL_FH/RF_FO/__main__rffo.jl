using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

include("../../build_graph.jl")
include("../__functions__.jl")
include("../__rffo__.jl")

env = Gurobi.Env()

nbr_thread = 8
silent = true 
graph_reduc = false
sp_time_limit = 60

# =============== Sur Nibi ===============
#instance = ARGS[1]
#graph_reduc = parse(Bool, ARGS[2])
# =============== En local sur mon ordi ===============
inst = ARGS[1]
instance = "../../INSTANCES/instances_literature_json/"*inst

# ============== Gestion des fichiers et la sortie =================
parts = split(instance, "/")
f_fold = parts[end-1]
Output_fold = "RESULTS_RF_FO/" * f_fold * "/"

# ============== Pré-run ================
file_essai = "../../INSTANCES/instances_json/A_MTN_2/127FL_5A_1.json"
inst_data = build_graph(file_essai)
path_f = "RESULTS_RF_FO/result_127FL_5A_1_RFFO.txt"
Otp_file = open(path_f, "w")
t_limit = 15
mdl = buildM(env, inst_data, nbr_thread, silent, graph_reduc, t_limit)
sz_day, ovlap = 3, 2
sol = relax_and_fix(mdl, inst_data, Otp_file, sz_day, ovlap)
sol_ = fix_and_optimize(mdl, sol.x, sol.y, inst_data, Otp_file, sz_day, ovlap)

# ============= Run officiel ============== 
inst_name = splitext(basename(instance))[1]

if graph_reduc
    path_file = Output_fold*"result_"*inst_name*"_RFFO_GR.txt"
else
    path_file = Output_fold*"result_"*inst_name*"_RFFO.txt"
end 

Output_file = open(path_file, "w")

write_both(Output_file, "======================= INSTANCE "*inst_name*"=======================")
write_both(Output_file, "\nGRAPH RÉDUCTION : "*string(graph_reduc)*"\n")

instance_data = build_graph(instance)
model = buildM(env, instance_data, nbr_thread, silent, graph_reduc, sp_time_limit)
size_day, overlap = 3, 2
write_both(Output_file, "RELAX AND FIX")
solution = relax_and_fix(model, instance_data, Output_file, size_day, overlap)
other_info = solution.other_info

rf_obj = solution.obj
rf_time = other_info["time"]
rf_nbr_iter = other_info["nbr_iter"]

write_both(Output_file, "\nRF OBJECTIVE VALUE : $rf_obj")
write_both(Output_file, "\nRF TIME  :  $rf_time seconds")

if solution.obj > 0
    write_both(Output_file, "FIX AND OPTIMIZE")
    sx = solution.x 
    sy = solution.y
    model = solution.other_info["model"]
    #= fo_sp_limit = 5
    set_optimizer_attribute(model, "TimeLimit", fo_sp_limit) 
    set_optimizer_attribute(model, "OutputFlag", 1) =#
    size_day, overlap = 3, 2
    solution = fix_and_optimize(model, sx, sy, instance_data, Output_file, size_day, overlap)
    other_info = solution.other_info

    fo_obj = solution.obj
    fo_time = solution.other_info["time"]
    fo_nbr_iter = other_info["nbr_iter"]

else
    write_both(Output_file, "NO FIX AND OPTIMIZE")
    fo_obj = solution.obj
    fo_time = 0.0
    fo_nbr_iter = 0
end 
rffo_time = rf_time + fo_time
rffo_nbr_iter = rf_nbr_iter + fo_nbr_iter
nbr_mtn = other_info["nbr_mtn"]
nbr_sts_used = other_info["nbr_stations_used"]

write_both(Output_file, "\nFO OBJECTIVE VALUE : $fo_obj")
write_both(Output_file, "\nFO TIME  : $fo_time seconds")
model = nothing
print_solution(solution, instance_data, Output_file, silent)

#close(Output_file)

dataframe = DataFrames.DataFrame(Instances = [inst_name], RF_Obj = [rf_obj], RF_Time = [rf_time], 
                                 RF_Iter = [rf_nbr_iter], FO_Obj = [fo_obj], FO_Time = [fo_time],
                                 FO_Iter = [fo_nbr_iter], RFFO_Time = [rffo_time], 
                                 RFFO_Iter = [rffo_nbr_iter], Nbr_mtn = [nbr_mtn], 
                                 Nbr_sts_used = nbr_sts_used)

if graph_reduc
    csv_summary_path = Output_fold * "summary_" * inst_name* "_RFFO_GR.csv"
else
    csv_summary_path = Output_fold * "summary_" * inst_name* "_RFFO.csv"
end 
# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

write_both(Output_file, "\nCSV Summary saved to: $csv_summary_path")
