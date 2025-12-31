using JuMP
using Gurobi
using JSON
using MathOptInterface
using DataFrames
using XLSX

include("model_amrp_new.jl")
include("benders_decomp_callback.jl")
include("functions.jl")
env = Gurobi.Env()
#= inst_name = ARGS[1]
preprocess = parse(Bool, ARGS[2])
=# 
instance = ARGS[1]

nbr_thread = 8
silent = true
time_limit = 7200
FH = true
FC = false
DY = false
MTN_CAP = false
#model_amrp(file*".json")

parts = split(instance, "/")
f_fold = parts[end-1]
Outputs_fold = "RESULTS/"*f_fold*"/"    

Instances, Obj, Obj_rho, Obj_lambda, Obj_phi, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], [], [], [], []
Time, Opt, Arc_reduc, Node_reduc, Nbr_mtn, Feasibilities = [], [], [], [], [], []
list_nbr_sts_used = []
nbr_feasibility_cuts, nbr_optimality_cuts, nbr_cuts, nbr_iterations = [], [], [], []
nbr_cover_cuts, MP_time, SP_time = [], [], []
inst_name = splitext(basename(instance))[1]
println(inst_name)
file_path = "INSTANCES/instances_json/A_MTN_2/"*inst_name
path_file = Outputs_fold*"result_"*inst_name*"_benders.txt"

Output_file = open(path_file, "w") 
write_both(Output_file, "-----------------------INSTANCE $inst_name--------------------------")
push!(Instances, inst_name)
solution = benders_decomp(env, instance, Output_file, nbr_thread, silent, time_limit)
push!(Obj, solution["obj"])
push!(MP_time, solution["mp_time"])
push!(SP_time, solution["sp_time"])
push!(Time, solution["time"])
push!(nbr_feasibility_cuts, solution["nbr_feasibility_cuts"])
push!(nbr_cover_cuts, solution["nbr_cover_cuts"])
push!(nbr_optimality_cuts, solution["nbr_optimality_cuts"])
push!(nbr_cuts, solution["nbr_cuts"])
push!(nbr_iterations, solution["nbr_iterations"])

write_both(Output_file, "\n=== Résultat final ===")
write_both(Output_file, "Objectif: $(solution["obj"])")
write_both(Output_file, "Temps total du master: $(round(solution["mp_time"], digits=2))s")
write_both(Output_file, "Temps total du sp: $(round(solution["sp_time"], digits=2))s") 
write_both(Output_file, "Temps total: $(solution["time"])s")
write_both(Output_file, "Nombre de coupes: $(solution["nbr_cuts"])")
write_both(Output_file, "Nombre de coupes faisabilité: $(solution["nbr_feasibility_cuts"])")
write_both(Output_file, "Nombre de cover inequalities: $(solution["nbr_cover_cuts"])")
write_both(Output_file, "Nombre de coupes d'optimalité: $(solution["nbr_optimality_cuts"])")
write_both(Output_file, "Nombre d'itérations: $(solution["nbr_iterations"])")
dataframe = DataFrames.DataFrame(Instances = Instances, UB = Obj, Time = Time, MP_Time = MP_time, SP_Time = SP_time,
                                Nbr_iterations = nbr_iterations, Nbr_feasibility_cuts = nbr_feasibility_cuts,
                                Nbr_optimality_cuts = nbr_optimality_cuts, Nbr_cover_cuts = nbr_cover_cuts, 
                                Nbr_cuts = nbr_cuts)
XLSX.writetable(Outputs_fold*"RESULTS_"*instance*"benders.xlsx", dataframe, overwrite=true)
