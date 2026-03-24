using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

#include("model_amrp.jl")
include("model_amrp_new.jl")
include("functions.jl")
include("../build_graph.jl")

env = Gurobi.Env()

instance = ARGS[1]
graph_reduc = parse(Bool, ARGS[2])
option = "MIN"

nbr_thread = 8
silent = false
time_limit = 43200
FH = true
FC = true
DY = true
MTN_CAP = false

# --- Gestion des fichiers de sortie ---
parts = split(instance, "/")
f_fold = parts[end-1]
Outputs_fold = "RESULTS_NOR/" * f_fold * "/"

inst_name = splitext(basename(instance))[1]

if graph_reduc
    path_file = Outputs_fold * "result_" * inst_name * "_GR.txt"
else
    path_file = Outputs_fold * "result_" * inst_name * ".txt"
end

Output_file = open(path_file, "w")
write_both(Output_file, "\nGRAPH REDUCING : " * string(graph_reduc))

write_both(Output_file, "-----------------------INSTANCE $inst_name--------------------------")
close(Output_file)

# --- Résolution du modèle ---
instance_data = build_graph(instance)
solution = model_amrp(env, instance_data, path_file, nbr_thread, silent, graph_reduc, time_limit, FH, FC, DY, MTN_CAP)
        
Output_file = open(path_file, "a")

other_info = solution.other_info
write_both(Output_file, "AVANT POSTPROCESSING")
write_both(Output_file, "obj_rho = " * string(other_info["obj_rho"]))
write_both(Output_file, "obj_lambda = " * string(other_info["obj_lambda"]))
write_both(Output_file, "obj_phi = " * string(other_info["obj_phi"]))
Feasibility = 1
#print_solution(solution, Output_file, silent)

write_both(Output_file, "APRÈS POSTPROCESSING")
solution = compute_indicators(solution, instance_data, FH, FC, DY)
write_both(Output_file, "obj_rho = " * string(other_info["obj_rho"]))
write_both(Output_file, "obj_lambda = " * string(other_info["obj_lambda"]))
write_both(Output_file, "obj_phi = " * string(other_info["obj_phi"]))
Feasibility = solution.other_info["feasible"]
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

dataframe = DataFrames.DataFrame(Instances = [inst_name], Rho = [Obj_rho], 
    Lambda = [Obj_lambda], Phi = [Obj_phi], UB = [Obj], LB = [Dual_obj], Gap = Gap, Nodes = Nbr_nodes, 
    Time = [Time], Feasible = [Feasibility], Opt = [Opt], Nbr_mtn = [Nbr_mtn], Nbr_sts_used = [Nbr_sts_used])


# --- Définir le chemin des fichiers CSV ---
if graph_reduc
    csv_summary_path = Outputs_fold * "summary_" * inst_name * "_MILP_GR.csv"
else
    csv_summary_path = Outputs_fold * "summary_" * inst_name * "_MILP.csv"
end


# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

write_both(Output_file, "CSV Summary saved to: "*string(csv_summary_path))
