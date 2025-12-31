using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

#include("model_amrp.jl")
include("../model_amrp_new.jl")
include("../functions.jl")

env = Gurobi.Env()

instance = ARGS[1]
preprocess = parse(Bool, ARGS[2])
option = "MIN"

nbr_thread = 8
silent = false
time_limit = 7200
FH = true
FC = false
DY = false
MTN_CAP = false

# --- Gestion des fichiers de sortie ---
parts = split(instance, "/")
f_fold = parts[end-1]
Output_fold = "RESULTS_NOR/" * f_fold * "/"

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
solution = model_amrp(env, instance, path_file, nbr_thread, silent, preprocess, time_limit, FH, FC, DY, MTN_CAP, option, option)

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

fh_parts = String[]
tk_parts = String[]
fd_parts = String[]

for ac in keys(mtn_infos)
    push!(fh_parts, "$ac : $(mtn_infos[ac][1])")
    push!(tk_parts, "$ac : $(mtn_infos[ac][2])")
    push!(fd_parts, "$ac : $(mtn_infos[ac][3])")
end

fh_ac = join(fh_parts, " | ")
tk_ac = join(tk_parts, " | ")
fd_ac = join(fd_parts, " | ")

#feasible = 1
Obj =  solution["obj"]
Time = solution["time"]
Obj_rho = solution["obj_rho"]
Obj_lambda = solution["obj_lambda"]
Obj_phi = solution["obj_phi"]
Dual_obj = solution["dual_obj"]
Gap = solution["gap"]
Nbr_nodes = solution["nbr_nodes"]
Arc_reduc = solution["arc_reduc"]
Node_reduc = solution["node_reduc"]
Nbr_mtn = solution["nbr_mtn"]
Nbr_sts_used = solution["nbr_sts_used"]
Feasible = solution["feasible"]
Fh_tk = solution["instance"]["fh_tk"]
Fh_day = solution["instance"]["fh_day"]

if solution["status"] == MOI.OPTIMAL
    Opt = 1
else
    Opt = 0
end

dataframe = DataFrames.DataFrame(Instances = [inst_name], fh_tk = [Fh_tk], fh_day = [Fh_day], Rho = [Obj_rho], 
    Lambda = [Obj_lambda], Phi = [Obj_phi], UB = [Obj], LB = [Dual_obj], Gap = Gap, Nodes = Nbr_nodes, 
    Time = [Time], Feasible = [Feasible], Opt = [Opt], Nbr_mtn = [Nbr_mtn], Nbr_sts_used = [Nbr_sts_used])
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
