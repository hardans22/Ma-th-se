using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob


include("relax_and_fix.jl")
include("../build_graph_copy.jl")
include("functions.jl")
env = Gurobi.Env()

#= inst_name = ARGS[1]
preprocess = parse(Bool, ARGS[2])
 =# 
instance = "354FL_8A_"
inst_list = ["136FL_5A_2"]

preprocess = false

nbr_thread = 8
silent = false
time_limit = 1800
FH = true
graph_reduc = true
#model_amrp(file*".json")
fold_1 = "RESULTS_NOR/A_MTN_2/"
Outputs_fold = fold_1

#Pré-run
file_essai = "../INSTANCES/instances_json/A_MTN_1/31FL_1A_1.json"
inst_data = build_graph(file_essai)
path_f = Outputs_fold*"result_31FL_1A_1_RF.txt"
mdl = buildM(env, inst_data, path_f, nbr_thread, silent, preprocess, time_limit, "RF")
sz_day, ovlap = 3, 2
sol = relax_and_fix(mdl, inst_data, sz_day, ovlap)
Instances, Obj = [], []
Time = []
list_nbr_sts_used = []
for v in 30:30    
    inst_name = instance*string(v)
    inst_path = "../INSTANCES/instances_literature_json/A_MTN_5/"*inst_name
    
    path_file = Outputs_fold*"result_"*inst_name*"_RF.txt"
    Output_file = open(path_file, "w") 
    write(Output_file, "INSTANCE "*inst_name)
    println("\nGRAPH REDUCTION : ", graph_reduc)
    write(Output_file, "\nGRAPH REDUCTION : "*string(graph_reduc))
    
    println("-----------------------INSTANCE $inst_name--------------------------")
    #write(Output_file, "\n-----------------------INSTANCE "*inst_name*"--------------------------")
    push!(Instances, inst_name)
    close(Output_file)

    instance_file = inst_path*".json"
    instance_data = build_graph(instance_file)
    model = buildM(env, instance_data, path_file, nbr_thread, silent, graph_reduc, time_limit, "RF")
    size_day, overlap = 3, 2
    solution = relax_and_fix(model, instance_data, size_day, overlap)
    push!(Obj, solution.obj)
    push!(Time, solution.other_info["time"])
    println("\nFINAL OBJECTIVE VALUE : ", solution.obj)
    println("\nTIME  : ", solution.other_info["time"], " seconds")
    mtn_infos = print_solution(solution, instance_data, Output_file, silent)
end 

dataframe = DataFrames.DataFrame(Instances = Instances, UB = Obj, Time = Time)

csv_summary_path = Outputs_fold * "summary_" * instance* "_RF.csv"


# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

println("\nCSV Summary saved to: ", csv_summary_path)
