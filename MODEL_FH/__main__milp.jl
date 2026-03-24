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
fold_1 = "RESULTS/A_MTN_3/"
Outputs_fold = fold_1

for graph_reduc in [false]
    Instances, Obj, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], []
    Time, Opt,Nbr_mtn = [], [], []
    list_nbr_sts_used = []
    fh_ac, tk_ac, fd_ac = [], [], []
    for v in 11:20
        inst_name = instance*string(v)
        inst_path = "../INSTANCES/instances_literature_json/A_MTN_3/"*inst_name
        if graph_reduc
            path_file = Outputs_fold * "result_" * inst_name * "_MILP_GR.txt"
        else
            path_file = Outputs_fold * "result_" * inst_name * "_MILP.txt"
        end

        Output_file = open(path_file, "w") 
        write_both(Output_file, "======================= INSTANCE "*inst_name*"=======================")
        write_both(Output_file, "\nGRAPH RÉDUCTION : "*string(graph_reduc)*"\n")

        push!(Instances, inst_name)
        close(Output_file)
        instance_file = inst_path*".json"
        instance_data = build_graph(instance_file)
        solution = model_amrp(env, instance_data, path_file, nbr_thread, silent, graph_reduc, time_limit)
        Output_file = open(path_file, "a") 
        other_info = solution.other_info
        print_solution(solution, instance_data, Output_file, silent)

        push!(Obj, solution.obj)
        push!(Time, other_info["time"])
        push!(Dual_obj, other_info["dual_obj"])
        push!(Gap, other_info["gap"])
        push!(Nbr_nodes, other_info["nbr_nodes"])
        #= push!(Arc_reduc, solution["arc_reduc"])
        push!(Node_reduc, solution["node_reduc"]) =#
        push!(Nbr_mtn, other_info["nbr_mtn"])
        push!(list_nbr_sts_used, other_info["nbr_sts_used"])

        if other_info["status"] == MOI.OPTIMAL
            push!(Opt, 1)
        else
            push!(Opt, 0)
        end
        #verification(solution, instance)
    end 

    dataframe = DataFrames.DataFrame(Instances = Instances, #= Arc_reduced = Arc_reduc, Node_reduced = Node_reduc, =# 
                                        UB = Obj, LB = Dual_obj, Gap = Gap, Nodes = Nbr_nodes, Time = Time, Opt = Opt, 
                                        Nbr_mtn = Nbr_mtn, Nbr_sts_used = list_nbr_sts_used)
    
    if graph_reduc
        csv_summary_path = Outputs_fold * "summary_"*instance*"MILP_GR.csv"
    else
        csv_summary_path = Outputs_fold * "summary_"*instance*"MILP.csv"
    end


    # --- Exporter les DataFrames en CSV ---
    CSV.write(csv_summary_path, dataframe)

    println("\nCSV Summary saved to: ", csv_summary_path)
end