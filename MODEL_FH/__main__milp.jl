using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob


include("model_amrp.jl")
include("../build_graph.jl")
include("functions.jl")
env = Gurobi.Env()

instance = "354FL_8A_"
inst_list = ["136FL_5A_2"]

nbr_thread = 8
silent = false
time_limit = 300
fold_1 = "RESULTS/A_MTN_1/"
Outputs_fold = fold_1

for graph_reduc in [false]
    Instances, Obj, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], []
    Time, Opt,Nbr_mtn = [], [], []
    list_nbr_sts_used = []
    fh_ac, tk_ac, fd_ac = [], [], []
    for v in 1:10   
        inst_name = instance*string(v)
        inst_path = "../INSTANCES/instances_literature_json/A_MTN_1/"*inst_name
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