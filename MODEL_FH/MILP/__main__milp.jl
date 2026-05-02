using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

#include("model_amrp.jl")
include("__model__milp.jl")
include("../../build_graph.jl")
include("../__functions__.jl")

env = Gurobi.Env()

nbr_thread = 8
silent = false
graph_reduc = false
time_limit = 18000


# =============== Sur Nibi ===============
#instance = ARGS[1]
#graph_reduc = parse(Bool, ARGS[2])

# =============== En local sur mon ordi ===============
inst = ARGS[1]
instance = "../../INSTANCES/instances_literature_json/"*inst

# ============== Gestion des fichiers et la sortie =================
parts = split(instance, "/")
f_fold = parts[end-1]
Outputs_fold = "RESULTS_MILP/" * f_fold * "/"
inst_name = splitext(basename(instance))[1]


if graph_reduc
    path_file = Outputs_fold * "result_" * inst_name * "_MILP_GR.txt"
else
    path_file = Outputs_fold * "result_" * inst_name * "_MILP.txt"
end

Output_file = open(path_file, "w")
write_both(Output_file, "======================= INSTANCE "*inst_name*" =======================")
write_both(Output_file, "\nGRAPH RÉDUCTION : "*string(graph_reduc)*"\n")
close(Output_file)

# --- Résolution du modèle ---
instance_data = build_graph(instance)
solution = model_amrp(env, instance_data, path_file, nbr_thread, silent, graph_reduc, time_limit)
        
Output_file = open(path_file, "a")

other_info = solution.other_info

print_solution(solution, instance_data, Output_file, silent)

#feasible = 1
Obj =  solution.obj
Time = other_info["time"]

Dual_obj = other_info["dual_obj"]
Gap = other_info["gap"]
Nbr_nodes = other_info["nbr_nodes"]

Nbr_mtn = other_info["nbr_mtn"]
Nbr_sts_used = other_info["nbr_sts_used"]

if other_info["status"] == MOI.OPTIMAL
    Opt = 1
else
    Opt = 0
end

dataframe = DataFrames.DataFrame(Instances = [inst_name], UB = [Obj], LB = [Dual_obj], Gap = Gap, Nodes = Nbr_nodes, 
    Time = [Time], Opt = [Opt], Nbr_mtn = [Nbr_mtn], Nbr_sts_used = [Nbr_sts_used])


# --- Définir le chemin des fichiers CSV ---
if graph_reduc
    csv_summary_path = Outputs_fold * "summary_" * inst_name * "_MILP_GR.csv"
else
    csv_summary_path = Outputs_fold * "summary_" * inst_name * "_MILP.csv"
end


# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)


write_both(Output_file, "CSV Summary saved to: $csv_summary_path")