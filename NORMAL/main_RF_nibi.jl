using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob


include("relax_and_fix.jl")
include("../build_graph.jl")
include("functions.jl")
env = Gurobi.Env()

inst_name = ARGS[1]
graph_reduc = parse(Bool, ARGS[2])

nbr_thread = 8
silent = false
time_limit = 1800
#model_amrp(file*".json")
# --- Gestion des fichiers de sortie ---
parts = split(instance, "/")
f_fold = parts[end-1]
Output_fold = "RESULTS_NOR/" * f_fold * "/"


#Pré-run
file_essai = "../INSTANCES/instances_json/A_MTN_2/127FL_5A_1.json"
inst_data = build_graph(file_essai)
path_f = Outputs_fold*"result_127FL_5A_1_RF.txt"
Otp_file = open(path_f, "w")
mdl = buildM(env, inst_data, nbr_thread, silent, graph_reduc, time_limit, "RF")
sz_day, ovlap = 3, 2
sol = relax_and_fix(mdl, inst_data, Otp_file, sz_day, ovlap)

#Run officiel 
inst_name = splitext(basename(instance))[1]
println(inst_name)

if graph_reduc
    path_file = Outputs_fold*"result_"*inst_name*"_RF_graph_reduc.txt"
else
    path_file = Outputs_fold*"result_"*inst_name*"_RF.txt"
end 

Output_file = open(path_file, "w")

write(Output_file, "INSTANCE "*inst_name)
println("\nGRAPH REDUCTION : ", graph_reduc)
write(Output_file, "\nGRAPH RÉDUCTION : "*string(graph_reduc)*"\n")

println("-----------------------INSTANCE $inst_name--------------------------")
#write(Output_file, "\n-----------------------INSTANCE "*inst_name*"--------------------------")

instance_data = build_graph(instance)
rf_sp_limit = 15
model = buildM(env, instance_data, nbr_thread, silent, graph_reduc, time_limit, "RF")
size_day, overlap = 3, 2
solution = relax_and_fix(model, instance_data, Output_file, size_day, overlap)
#Output_file = open(path_file, "a")
Obj = solution.obj
Time = solution.other_info["time"]
println("\nFINAL OBJECTIVE VALUE : ", solution.obj)
println("\nTIME  : ", solution.other_info["time"], " seconds")
mtn_infos = print_solution(solution, instance_data, Output_file, silent)
model = nothing
#close(Output_file)


dataframe = DataFrames.DataFrame(Instances = [inst_name], UB = [Obj], Time = [Time])

if graph_reduc
    csv_summary_path = Outputs_fold * "summary_" * instance* "_RF_graph_reduc.csv"
else
    csv_summary_path = Outputs_fold * "summary_" * instance* "_RF.csv"
end 
# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

println("\nCSV Summary saved to: ", csv_summary_path)
