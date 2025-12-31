using JuMP, Gurobi, JSON, MathOptInterface, DataFrames, XLSX, CSV

include("../build_graph.jl")
include("benders_decomp_callback.jl")
include("../functions.jl")

env = Gurobi.Env()
#= inst_name = ARGS[1]
preproces
s = parse(Bool, ARGS[2])
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
Outputs_fold = "RESULTS_BENDERS/"*f_fold*"/"    

inst_name = splitext(basename(instance))[1]
println(inst_name)
path_file = Outputs_fold*"result_"*inst_name*"_benders.txt"

Output_file = open(path_file, "w") 
write_both(Output_file, "-----------------------INSTANCE $inst_name--------------------------")
Instances = inst_name
close(Output_file)

solution = benders_decomp(env, instance, path_file, nbr_thread, silent, time_limit)

Output_file = open(path_file, "a") 

Obj = solution["obj"]
MP_time = solution["mp_time"]
SP_time = solution["sp_time"]
prep1_time = solution["prep1_time"]
prep2_time = solution["prep2_time"]
Time = solution["time"]
nbr_feasibility_cuts = solution["nbr_feasibility_cuts"]
nbr_cover_cuts = solution["nbr_cover_cuts"]
nbr_optimality_cuts = solution["nbr_optimality_cuts"]
nbr_cuts = solution["nbr_cuts"]
nbr_iterations = solution["nbr_iterations"]

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

dataframe = DataFrames.DataFrame(Instances = [Instances], UB = [Obj], MP_Time = [MP_time], SP_Time = [SP_time],
                                Prep1_time = [prep1_time], Prep2_time = [prep2_time], Time = [Time], 
                                Nbr_iterations = [nbr_iterations], Nbr_feasibility_cuts = [nbr_feasibility_cuts],
                                Nbr_optimality_cuts = [nbr_optimality_cuts], Nbr_cover_cuts = [nbr_cover_cuts], 
                                Nbr_cuts = [nbr_cuts])

csv_summary_path = Outputs_fold*"summary_"*inst_name*"_benders.csv"
# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

println("CSV Summary saved to: ", csv_summary_path)
