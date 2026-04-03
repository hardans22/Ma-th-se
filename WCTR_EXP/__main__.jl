using JuMP, Gurobi, JSON, CSV, DataFrames, XLSX, Glob


include("model_amrp.jl")
include("../build_graph.jl")
include("__functions__.jl")
env = Gurobi.Env()

#= inst_name = ARGS[1]
preprocess = parse(Bool, ARGS[2])
 =# 
instance = "354FL_8A_"
preprocess = false

nbr_thread = 8
silent = false
time_limit = 20
obj_min = "FH"
#model_amrp(file*".json")
fold_1 = "RESULTS/A_MTN_5/"
Outputs_fold = fold_1

for graph_reduc in [false]
    Instances, Obj, Obj_rho, Obj_lambda, Obj_phi, Dual_obj = [], [], [], [], [], []
    Gap, Nbr_nodes, Time, Opt, Nbr_mtn = [], [], [], [], []
    list_nbr_sts_used = []
    for v in 21:30
        inst_name = instance*string(v)
        inst_path = "../INSTANCES/instances_literature_json/A_MTN_5/"*inst_name
        if graph_reduc
            path_file = Outputs_fold*"result_"*inst_name*"_graph_reduc.txt"
        else
            path_file = Outputs_fold*"result_"*inst_name*".txt"
        end
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
        solution = model_amrp(env, instance_data, path_file, nbr_thread, time_limit, obj_min, silent, graph_reduc)
        Output_file = open(path_file, "a") 
        other_info = solution.other_info
        #print_solution(solution, instance_data, Output_file, silent)
        
        push!(Obj, solution.obj)
        push!(Time, other_info["time"])
        push!(Obj_rho, other_info["obj_rho"])
        push!(Obj_lambda, other_info["obj_lambda"])
        push!(Obj_phi, other_info["obj_phi"])
        push!(Dual_obj, other_info["dual_obj"])
        push!(Gap, other_info["gap"])
        push!(Nbr_nodes, other_info["nbr_nodes"])
        push!(Nbr_mtn, other_info["nbr_mtn"])
        push!(list_nbr_sts_used, other_info["nbr_sts_used"])

        if other_info["status"] == MOI.OPTIMAL
            push!(Opt, 1)
        else
            push!(Opt, 0)
        end
        #verification(solution, instance)
    end 

    dataframe = DataFrames.DataFrame(Instances = Instances, Rho = Obj_rho, Lambda = Obj_lambda, Phi = Obj_phi,
                                        UB = Obj, LB = Dual_obj, Gap = Gap, Nodes = Nbr_nodes, Time = Time,
                                        Opt = Opt, Nbr_mtn = Nbr_mtn, Nbr_sts_used = list_nbr_sts_used)
    
    if graph_reduc
        csv_summary_path = Outputs_fold * "summary_" * instance * "_graph_reduc.csv"
    else
        csv_summary_path = Outputs_fold * "summary_" * instance * ".csv"
    end

    # --- Exporter les DataFrames en CSV ---
    CSV.write(csv_summary_path, dataframe)
    #CSV.write(csv_aircrafts_path, dataframe2)

    println("\nCSV Summary saved to: ", csv_summary_path)
    #println("CSV Aircrafts saved to: ", csv_aircrafts_path)
end