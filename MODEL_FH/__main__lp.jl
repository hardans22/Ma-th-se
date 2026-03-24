using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob


include("model_amrp.jl")
include("../build_graph.jl")
include("functions.jl")
env = Gurobi.Env()

#= inst_name = ARGS[1]
preprocess = parse(Bool, ARGS[2])
 =# 
instance = "354FL_8A_"
inst_list = ["136FL_5A_2"]
#= inst_list = [
    "354FL_8A_1", "354FL_8A_2", "354FL_8A_3", "354FL_8A_4", "354FL_8A_5",
    "354FL_8A_6", "354FL_8A_7", "354FL_8A_8", "354FL_8A_9", "354FL_8A_10",
    "354FL_8A_11", "354FL_8A_12", "354FL_8A_13", "354FL_8A_14", "354FL_8A_15",
    "354FL_8A_16", "354FL_8A_17", "354FL_8A_18", "354FL_8A_19", "354FL_8A_20",
    "354FL_8A_21", "354FL_8A_22", "354FL_8A_23", "354FL_8A_24", "354FL_8A_25",
    "354FL_8A_26", "354FL_8A_27", "354FL_8A_28", "354FL_8A_29", "354FL_8A_30"
]
 =#
preprocess = false

nbr_thread = 8
silent = false
time_limit = 300
fold_1 = "RESULTS/A_MTN_1/"
Outputs_fold = fold_1

for graph_reduc in [false]
    Instances, Obj, Time = [], [], []
    for v in 1:10
        inst_name = instance*string(v)
        inst_path = "../INSTANCES/instances_literature_json/A_MTN_1/"*inst_name
        if graph_reduc
            path_file = Outputs_fold * "result_" * inst_name * "_LP_GR.txt"
        else
            path_file = Outputs_fold * "result_" * inst_name * "_LP.txt"
        end

        Output_file = open(path_file, "w") 
        write_both(Output_file, "======================= INSTANCE "*inst_name*"=======================")
        write_both(Output_file, "\nGRAPH RÉDUCTION : "*string(graph_reduc)*"\n")
        push!(Instances, inst_name)
        
        instance_file = inst_path*".json"
        instance_data = build_graph(instance_file)
        solution = model_amrp_lp(env, instance_data, nbr_thread, silent, graph_reduc)
        other_info = solution.other_info

        push!(Obj, solution.obj)
        push!(Time, other_info["time"])
    end 

    dataframe = DataFrames.DataFrame(Instances = Instances, Obj = Obj, Time = Time)
    
    if graph_reduc
        csv_summary_path = Outputs_fold * "summary_"*instance*"LP_GR.csv"
    else
        csv_summary_path = Outputs_fold * "summary_"*instance*"LP.csv"
    end


    # --- Exporter les DataFrames en CSV ---
    CSV.write(csv_summary_path, dataframe)

    println("\nCSV Summary saved to: ", csv_summary_path)
end