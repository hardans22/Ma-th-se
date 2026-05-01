using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

#include("model_amrp.jl")
include("model_amrp.jl")
include("functions.jl")
include("../build_graph.jl")

env = Gurobi.Env()

instance = ARGS[1]
graph_reduc = parse(Bool, ARGS[2])
option = "MIN"

nbr_thread = 8
silent = false
time_limit = 120

parts = split(instance, "/")
f_fold = parts[end-1]
Outputs_fold = "RESULTS/" * f_fold * "/"

inst_name = splitext(basename(instance))[1]

if graph_reduc
    path_file = Outputs_fold * "result_" * inst_name * "_LP_GR.txt"
else
    path_file = Outputs_fold * "result_" * inst_name * "_LP.txt"
end

Output_file = open(path_file, "w")
write_both(Output_file, "======================= INSTANCE "*inst_name*"=======================")
write_both(Output_file, "\nGRAPH RÉDUCTION : "*string(graph_reduc)*"\n")

# --- Résolution du modèle ---
instance_data = build_graph(instance)
solution = model_amrp_lp(env, instance_data, nbr_thread, silent, graph_reduc)
other_info = solution.other_info
Obj =  solution.obj
Time = other_info["time"]

dataframe = DataFrames.DataFrame(Instances = [inst_name], UB = [Obj], Time = [Time])

# --- Définir le chemin des fichiers CSV ---
if graph_reduc
    csv_summary_path = Outputs_fold * "summary_" * inst_name * "_LP_GR.csv"
else
    csv_summary_path = Outputs_fold * "summary_" * inst_name * "_LP.csv"
end

# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

write_both(Output_file, "CSV Summary saved to: $csv_summary_path")