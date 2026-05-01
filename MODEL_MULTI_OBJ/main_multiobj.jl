using JuMP, Gurobi, JSON, CSV, DataFrames, XLSX, Glob


include("amrp_multiobj.jl")
include("../build_graph.jl")
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
time_limit = 600
FH = true
FC = false
DY = false
option = "MEAN"
MTN_CAP = true
#model_amrp(file*".json")
fold_1 = "RESULTS_NOR/A_MTN_1/"
Outputs_fold = fold_1

for graph_reduc in [false]
    Instances, Obj, Obj_rho, Obj_lambda, Obj_phi, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], [], [], [], []
    Time, Opt, Arc_reduc, Node_reduc, Nbr_mtn, Feasibilities = [], [], [], [], [], []
    list_nbr_sts_used = []
    fh_ac, tk_ac, fd_ac = [], [], []
    for v in 1:1
        inst_name = instance*string(v)
        inst_path = "../INSTANCES/instances_literature_json/A_MTN_1/"*inst_name
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
        solution = model_amrp(env, instance_data, path_file, nbr_thread, silent, graph_reduc, time_limit)
        Output_file = open(path_file, "a") 
        other_info = solution.other_info
        print_solution(solution, instance_data, Output_file, silent)

        println("AVANT POSTPROCESSING")
        println("obj_rho = ", other_info["obj_rho"])
        println("obj_lambda = ", other_info["obj_lambda"])
        println("obj_phi = ", other_info["obj_phi"])
        other_info["feasible"] = 1
        
        push!(Obj, solution.obj)
        push!(Time, other_info["time"])
        push!(Obj_rho, other_info["obj_rho"])
        push!(Obj_lambda, other_info["obj_lambda"])
        push!(Obj_phi, other_info["obj_phi"])
        push!(Dual_obj, other_info["dual_obj"])
        push!(Gap, other_info["gap"])
        push!(Nbr_nodes, other_info["nbr_nodes"])
        #= push!(Arc_reduc, solution["arc_reduc"])
        push!(Node_reduc, solution["node_reduc"]) =#
        push!(Nbr_mtn, other_info["nbr_mtn"])
        push!(list_nbr_sts_used, other_info["nbr_sts_used"])
        push!(Feasibilities, other_info["feasible"])

        if other_info["status"] == MOI.OPTIMAL
            push!(Opt, 1)
        else
            push!(Opt, 0)
        end
        #verification(solution, instance)
    end 

    dataframe = DataFrames.DataFrame(Instances = Instances, #= Arc_reduced = Arc_reduc, Node_reduced = Node_reduc, =# 
                                        Rho = Obj_rho, Lambda = Obj_lambda, Phi = Obj_phi, UB = Obj, LB = Dual_obj, 
                                        Gap = Gap, Nodes = Nbr_nodes, Time = Time, Feasible = Feasibilities, Opt = Opt, 
                                        Nbr_mtn = Nbr_mtn, Nbr_sts_used = list_nbr_sts_used)
    #dataframe2 = DataFrames.DataFrame(Instances = Instances, FH = fh_ac, FC = tk_ac, DY = fd_ac)
    
    if graph_reduc
        csv_summary_path = Outputs_fold * "summary_" * instance * "_graph_reduc.csv"
        csv_aircrafts_path = Outputs_fold * "aircrafts_" * instance * "_graph_reduc.csv"
    else
        csv_summary_path = Outputs_fold * "summary_" * instance * ".csv"
        csv_aircrafts_path = Outputs_fold * "aircrafts_" * instance * ".csv"
    end

    # --- Exporter les DataFrames en CSV ---
    CSV.write(csv_summary_path, dataframe)
    #CSV.write(csv_aircrafts_path, dataframe2)

    println("\nCSV Summary saved to: ", csv_summary_path)
    #println("CSV Aircrafts saved to: ", csv_aircrafts_path)
end