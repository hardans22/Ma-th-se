using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob


include("relax_and_fix.jl")
include("../build_graph.jl")
include("functions.jl")
env = Gurobi.Env()

instance = ARGS[1]
graph_reduc = parse(Bool, ARGS[2])

nbr_thread = 8
silent = false
#model_amrp(file*".json")
# --- Gestion des fichiers de sortie ---
parts = split(instance, "/")
f_fold = parts[end-1]
Output_fold = "RESULTS_NOR/" * f_fold * "/"

#Pré-run
file_essai = "../INSTANCES/instances_json/A_MTN_2/127FL_5A_1.json"
inst_data = build_graph(file_essai)
path_f = "RESULTS_NOR/A_MTN_2/result_127FL_5A_1_RF.txt"
Otp_file = open(path_f, "w")
t_limit = 15
mdl = buildM(env, inst_data, nbr_thread, silent, graph_reduc, t_limit, "RF")
sz_day, ovlap = 3, 2
sol = relax_and_fix(mdl, inst_data, Otp_file, sz_day, ovlap)

#Run officiel 
inst_name = splitext(basename(instance))[1]
println(inst_name)

if graph_reduc
    path_file = Output_fold*"result_"*inst_name*"_RF_graph_reduc.txt"
else
    path_file = Output_fold*"result_"*inst_name*"_RF.txt"
end 

Output_file = open(path_file, "w")

write(Output_file, "INSTANCE "*inst_name)
println("\nGRAPH REDUCTION : ", graph_reduc)
write(Output_file, "\nGRAPH RÉDUCTION : "*string(graph_reduc)*"\n")

println("-----------------------INSTANCE $inst_name--------------------------")
#write(Output_file, "\n-----------------------INSTANCE "*inst_name*"--------------------------")

instance_data = build_graph(instance)
time_limit = 60
model = buildM(env, instance_data, nbr_thread, silent, graph_reduc, time_limit, "RF")
size_day, overlap = 3, 2
println("RELAX AND FIX")
solution = relax_and_fix(model, instance_data, Output_file, size_day, overlap)
rf_obj = solution.obj
rf_time = solution.other_info["time"]
println("\nFINAL OBJECTIVE VALUE : ", rf_obj)
println("\nTIME  : ", rf_time, " seconds")

if solution.obj > 0
    println("FIX AND OPTIMIZE")
    sx = solution.x 
    model = solution.other_info["model"]
    fo_sp_limit = 5
    #= set_optimizer_attribute(model, "TimeLimit", fo_sp_limit) 
    set_optimizer_attribute(model, "OutputFlag", 1) =#
    size_day, overlap = 3, 2
    solution = fix_and_optimize(model, sx, instance_data, Output_file, size_day, overlap)

    fo_obj = solution.obj
    fo_time = solution.other_info["time"]
else
    println("NO FIX AND OPTIMIZE")
    fo_obj = solution.obj
    fo_time = 0.0
end 
println("\nFINAL OBJECTIVE VALUE : ", fo_obj)
println("\nTIME  : ", fo_time, " seconds")
model = nothing
print_solution(solution, instance_data, Output_file, silent)

#close(Output_file)


dataframe = DataFrames.DataFrame(Instances = [inst_name], RF_Obj = [rf_obj], RF_Time = [rf_time], FO_Obj = [fo_obj], FO_Time = [fo_time])

if graph_reduc
    csv_summary_path = Output_fold * "summary_" * inst_name* "_RF_graph_reduc.csv"
else
    csv_summary_path = Output_fold * "summary_" * inst_name* "_RF.csv"
end 
# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

println("\nCSV Summary saved to: ", csv_summary_path)
