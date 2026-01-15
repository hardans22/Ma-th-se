using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob


include("model_amrp_new.jl")
include("../build_graph_copy.jl")
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
time_limit = 1800
FH = true
FC = false
DY = false
option = "MEAN"
MTN_CAP = false
#model_amrp(file*".json")
fold_1 = "RESULTS_NOR/A_MTN_2/"
Outputs_fold = fold_1


for graph_reduc in [true]
    Instances, Obj, Obj_rho, Obj_lambda, Obj_phi, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], [], [], [], []
    Time, Opt, Arc_reduc, Node_reduc, Nbr_mtn, Feasibilities = [], [], [], [], [], []
    list_nbr_sts_used = []
    fh_ac, tk_ac, fd_ac = [], [], []
    for v in 30:30
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
        solution = model_amrp(env, instance_data, path_file, nbr_thread, silent, graph_reduc, time_limit, FH, FC, DY, MTN_CAP, option, option)
        Output_file = open(path_file, "a") 
        other_info = solution.other_info
        println("AVANT POSTPROCESSING")
        println("obj_rho = ", other_info["obj_rho"])
        println("obj_lambda = ", other_info["obj_lambda"])
        println("obj_phi = ", other_info["obj_phi"])
        other_info["feasible"] = 1
        
        
        #= println("APRÈS POSTPROCESSING")
        solution = compute_indicators(solution, FH, FC, DY)
        println("obj_rho = ", solution["obj_rho"])
        println("obj_lambda = ", solution["obj_lambda"])
        println("obj_phi = ", solution["obj_phi"])
         =#
        mtn_infos = print_solution(solution, instance_data, Output_file, silent)
        
        #= fh_ac_str, tk_ac_str, fd_ac_str = "", "", ""
        for ac in keys(mtn_infos)
            fh_ac_str *= ac*" : "*string(mtn_infos[ac][1])*" | "
            tk_ac_str *= ac*" : "*string(mtn_infos[ac][2])*" | "     
            fd_ac_str *= ac*" : "*string(mtn_infos[ac][3])*" | "     
        end
        push!(fh_ac, fh_ac_str)
        push!(tk_ac, tk_ac_str)
        push!(fd_ac, fd_ac_str)
         =#
        #solution["feasible"] = 1
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