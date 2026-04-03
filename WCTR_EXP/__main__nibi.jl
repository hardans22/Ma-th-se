using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

#include("model_amrp.jl")
include("model_amrp.jl")
include("__functions__.jl")
include("../build_graph.jl")

env = Gurobi.Env()

instance = ARGS[1]
graph_reduc = parse(Bool, ARGS[2])

nbr_thread = 8
silent = false
time_limit = 18000
obj_min = "FH"

# --- Gestion des fichiers de sortie ---
parts = split(instance, "/")
f_fold = parts[end-1]
Outputs_fold = "RESULTS/" * f_fold * "/"

inst_name = splitext(basename(instance))[1]

if graph_reduc
    path_file = Outputs_fold * "result_" * inst_name * "_MILP_GR.txt"
else
    path_file = Outputs_fold * "result_" * inst_name * "_MILP.txt"
end

Output_file = open(path_file, "w")
write_both(Output_file, "\nGRAPH REDUCING : " * string(graph_reduc))

write_both(Output_file, "-----------------------INSTANCE $inst_name--------------------------")
close(Output_file)

# --- Résolution du modèle ---
instance_data = build_graph(instance)
solution = model_amrp(env, instance_data, path_file, nbr_thread, time_limit, obj_min, silent, graph_reduc)
        
Output_file = open(path_file, "a")

other_info = solution.other_info

write_both(Output_file, "Remaining value of each indicator")
write_both(Output_file, "obj_rho = " * string(other_info["obj_rho"]))
write_both(Output_file, "obj_lambda = " * string(other_info["obj_lambda"]))
write_both(Output_file, "obj_phi = " * string(other_info["obj_phi"]))
print_solution(solution, instance_data, Output_file, silent)

Obj =  solution.obj
Time = other_info["time"]

Dual_obj = other_info["dual_obj"]
Gap = other_info["gap"]
Nbr_nodes = other_info["nbr_nodes"]

Nbr_mtn = other_info["nbr_mtn"]
Nbr_sts_used = other_info["nbr_sts_used"]
Obj_rho = other_info["obj_rho"]
Obj_lambda = other_info["obj_lambda"]
Obj_phi = other_info["obj_phi"]
        
if other_info["status"] == MOI.OPTIMAL
    Opt = 1
else
    Opt = 0
end

dataframe = DataFrames.DataFrame(Instances = [inst_name], Rho = [Obj_rho], Lambda = [Obj_lambda],
                    Phi = [Obj_phi], UB = [Obj], LB = [Dual_obj], Gap = Gap, Nodes = Nbr_nodes, 
                    Time = [Time], Opt = [Opt], Nbr_mtn = [Nbr_mtn], Nbr_sts_used = [Nbr_sts_used])


# --- Définir le chemin des fichiers CSV ---
if graph_reduc
    csv_summary_path = Outputs_fold * "summary_" * inst_name * "_MILP_GR.csv"
else
    csv_summary_path = Outputs_fold * "summary_" * inst_name * "_MILP.csv"
end


# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

write_both(Output_file, "CSV Summary saved to: "*string(csv_summary_path))
