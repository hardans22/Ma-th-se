using JuMP
using Gurobi
using JSON
using MathOptInterface
using DataFrames
using XLSX

include("model_amrp.jl")
#include("base_model_amrp.jl")
include("functions.jl")
preprocess = parse(Bool, ARGS[1])

nbr_thread = 10
silent = false
time_limit = 14400
#model_amrp(file*".json")

Outputs_fold = "RESULTS/"
for preprocess in [false, true]
    list_instances, list_obj, list_dual_obj, list_gap, list_nbr_nodes  = [], [], [], [], []
    list_time, list_opt, list_arc_reduc, list_node_reduc, list_nbr_mtn = [], [], [], [], []
    list_nbr_sts_used = []
    println("\nPREPROCESSING : ", preprocess)
    for inst_name in ["75FL_7A_3D", "76FL_7A_3D", "78FL_7A_3D", "86FL_7A_3D", "89FL_7A_3D"]
        file_path = "instances_json/"*inst_name
        println("INSTANCE : $inst_name")
        if preprocess
            path_file = Outputs_fold*"result_"*inst_name*"_with_prep.txt"
        else 
            path_file = Outputs_fold*"result_"*inst_name*"_without_prep.txt"
        end
        Output_file = open(path_file, "w") 
        write(Output_file, "INSTANCE "*inst_name)
        write(Output_file, "\nPREPROCESSING : "*string(preprocess))
        for i in 1:10
            println("-----------------------INSTANCE $i--------------------------")
            write(Output_file, "\n-----------------------INSTANCE "*string(i)*"--------------------------")
            push!(list_instances, inst_name*"_"*string(i))
            solution = model_amrp(file_path*"_"*string(i)*".json", nbr_thread, silent, preprocess, time_limit)
            push!(list_obj, solution["obj"])
            push!(list_dual_obj, solution["dual_obj"])
            push!(list_gap, solution["gap"])
            push!(list_nbr_nodes, solution["nbr_nodes"])
            push!(list_time, solution["time"])
            push!(list_arc_reduc, solution["arc_reduc"])
            push!(list_node_reduc, solution["node_reduc"])
            push!(list_nbr_mtn, solution["nbr_mtn"])
            push!(list_nbr_sts_used, solution["nbr_sts_used"])
            if solution["status"] == MOI.OPTIMAL
                push!(list_opt, 1)
            else
                push!(list_opt, 0)
            end
            print_solution(solution, Output_file, silent)
        end 
    end
    dataframe = DataFrames.DataFrame(Instances = list_instances, Arc_reduced = list_arc_reduc, 
                                    Node_reduced = list_node_reduc, UB = list_obj, LB = list_dual_obj, 
                                    Gap = list_gap, Nodes = list_nbr_nodes, Time = list_time, Opt = list_opt, 
                                    Nbr_mtn = list_nbr_mtn, Nbr_sts_used = list_nbr_sts_used)
    if preprocess
        XLSX.writetable(Outputs_fold*"result_with_prep.xlsx", dataframe, overwrite=true)
    else
        XLSX.writetable(Outputs_fold*"result_without_prep.xlsx", dataframe, overwrite=true)
    end
end