using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

#include("model_amrp.jl")
include("model_amrp_new.jl")
include("functions.jl")

env = Gurobi.Env()

instance = ARGS[1]
preprocess = parse(Bool, ARGS[2])
option = ARGS[3]
#= instance = "./instances_json/104FL_10A_3D_5MS_1.json"
preprocess = false
 =#
nbr_thread = 8
silent = true
time_limit = 7200
FH = true
FC = true
DY = false
MTN_CAP = false

# --- Gestion des fichiers de sortie ---
parts = split(instance, "/")
f_fold = parts[end-1]
fold_1 = "RESULTS/" * f_fold * "/"

fold_2 = ""
if FH && FC && option == "MIN" 
    fold_2 = "FH_FC_MIN/"
elseif FH && FC && option == "MAX" 
    fold_2 = "FH_FC_MAX/"
elseif FH && FC && option == "MEDIAN" 
    fold_2 = "FH_FC_MEDIAN/"
elseif FH && FC && option == "MEAN" 
    fold_2 = "FH_FC_MEAN/"
elseif FH && DY && option == "MIN" 
    fold_2 = "FH_DY_MIN/"
elseif FH && DY && option == "MAX" 
    fold_2 = "FH_DY_MAX/"
elseif FH && DY && option == "MEDIAN" 
    fold_2 = "FH_DY_MEDIAN/"
elseif FH && DY && option == "MEAN" 
    fold_2 = "FH_DY_MEAN/"
end

Outputs_fold = fold_1*fold_2


inst_name = splitext(basename(instance))[1]
println(inst_name)

if preprocess
    path_file = Outputs_fold * "result_" * inst_name * "_with_prep.txt"
else
    path_file = Outputs_fold * "result_" * inst_name * "_without_prep.txt"
end

Output_file = open(path_file, "w")
write(Output_file, "INSTANCE "*inst_name)
println("\nPREPROCESSING : ", preprocess)
write(Output_file, "\nPREPROCESSING : "*string(preprocess))

println("-----------------------INSTANCE $inst_name--------------------------")
write(Output_file, "\n-----------------------INSTANCE "*string(inst_name)*"--------------------------")
close(Output_file)

# --- Résolution du modèle ---
solution = model_amrp(env, instance, Output_file, nbr_thread, silent, preprocess, time_limit, FH, FC, DY, MTN_CAP, option,option)

Output_file = open(path_file, "a")

println("AVANT CALCUL DES AUTRES INDICATEURS")
println("obj_rho = ", solution["obj_rho"])
println("obj_lambda = ", solution["obj_lambda"])
println("obj_phi = ", solution["obj_phi"])
write(Output_file, "\nAVANT CALCUL DES AUTRES INDICATEURS")
write(Output_file, "\nObj_rho = "*string(solution["obj_rho"]))
write(Output_file, "\nObj_lambda = "*string(solution["obj_lambda"]))
write(Output_file, "\nObj_phi = "*string(solution["obj_phi"]))

#feasible = 1    
solution = compute_indicators(solution, FH, FC, DY)
println("APRÈS CALCUL DES AUTRES INDICATEURS")
println("obj_rho = ", solution["obj_rho"])
println("obj_lambda = ", solution["obj_lambda"])
println("obj_phi = ", solution["obj_phi"])
write(Output_file, "\nAVANT CALCUL DES AUTRES INDICATEURS")
write(Output_file, "\nObj_rho = "*string(solution["obj_rho"]))
write(Output_file, "\nObj_lambda = "*string(solution["obj_lambda"]))
write(Output_file, "\nObj_phi = "*string(solution["obj_phi"]))

mtn_infos = print_solution(solution, Output_file, silent)
fh_ac, tk_ac, fd_ac = "", "", ""
for ac in keys(mtn_infos)
    fh_ac *= ac*" : "*string(mtn_infos[ac][1])*" | "
    tk_ac *= ac*" : "*string(mtn_infos[ac][2])*" | "     
    fd_ac *= ac*" : "*string(mtn_infos[ac][3])*" | "     
end

obj = solution["obj"]
obj_rho = solution["obj_rho"]
obj_lambda = solution["obj_lambda"]
obj_phi = solution["obj_phi"]
dual_obj = solution["dual_obj"]
gap = solution["gap"]
nbr_nodes = solution["nbr_nodes"]
time = solution["time"]
arc_reduc = solution["arc_reduc"]
node_reduc = solution["node_reduc"]
nbr_mtn = solution["nbr_mtn"]
nbr_sts_used = solution["nbr_sts_used"]
opt = solution["status"] == MOI.OPTIMAL ? 1 : 0
feasible = solution["feasible"]
fh_tk = instance["fh_tk"]
fh_day = instance["fh_day"]

dataframe = DataFrames.DataFrame( Instances = [inst_name], fh_tk = [fh_tk], fh_day = [fh_day], #= Arc_reduced = [arc_reduc], Node_reduced = [node_reduc], =# Rho = [obj_rho],
    Lambda = [obj_lambda], Phi = [obj_phi], UB = [obj], LB = [dual_obj], Gap = [gap], Nodes = [nbr_nodes],
    Time = [time], Opt = [opt], Feasible = [feasible], Nbr_mtn = [nbr_mtn], Nbr_sts_used = [nbr_sts_used]
)

dataframe2 = DataFrames.DataFrame(Instances = [inst_name], FH = [fh_ac], FC = [tk_ac], DY = [fd_ac])

# --- Définir le chemin des fichiers CSV ---
if preprocess
    csv_summary_path = Outputs_fold * "summary_" * inst_name * "_with_prep.csv"
    csv_aircrafts_path = Outputs_fold * "aircrafts_" * inst_name * "_with_prep.csv"
else
    csv_summary_path = Outputs_fold * "summary_" * inst_name * "_without_prep.csv"
    csv_aircrafts_path = Outputs_fold * "aircrafts_" * inst_name * "_without_prep.csv"
end

# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)
CSV.write(csv_aircrafts_path, dataframe2)

println("CSV Summary saved to: ", csv_summary_path)
println("CSV Aircrafts saved to: ", csv_aircrafts_path)

