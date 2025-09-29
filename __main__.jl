using JuMP
using Gurobi
using JSON
using MathOptInterface
using DataFrames
using XLSX

include("model_amrp.jl")
#include("base_model_amrp.jl")
include("functions.jl")
#= inst_name = ARGS[1]
preprocess = parse(Bool, ARGS[2])
 =# 
inst_name = "358FL_15A_1"
#preprocess = false

file_path = "instances_json/A_MTN_3/"*inst_name
nbr_thread = 10
silent = false
time_limit = 20
#model_amrp(file*".json")

Outputs_fold = "RESULTS/"
println("INSTANCE : $inst_name")

Instances, Obj, Obj_rho, Obj_lambda, Obj_phi, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], [], [], [], []
Time, Opt, Arc_reduc, Node_reduc, Nbr_mtn = [], [], [], [], []
list_nbr_sts_used = []
for preprocess in [false]
    if preprocess
        path_file = Outputs_fold*"result_"*inst_name*"_with_prep.txt"
    else 
        path_file = Outputs_fold*"result_"*inst_name*"_without_prep.txt"
    end
    Output_file = open(path_file, "w") 
    write(Output_file, "INSTANCE "*inst_name)
    println("\nPREPROCESSING : ", preprocess)
    write(Output_file, "\nPREPROCESSING : "*string(preprocess))
    
    println("-----------------------INSTANCE $inst_name--------------------------")
    write(Output_file, "\n-----------------------INSTANCE "*inst_name*"--------------------------")
    push!(Instances, inst_name)
    solution = model_amrp(file_path*".json", nbr_thread, silent, preprocess, time_limit)
    push!(Obj, solution["obj"])
    push!(Obj_rho, solution["obj_rho"])
    push!(Obj_lambda, solution["obj_lambda"])
    push!(Obj_phi, solution["obj_phi"])
    push!(Dual_obj, solution["dual_obj"])
    push!(Gap, solution["gap"])
    push!(Nbr_nodes, solution["nbr_nodes"])
    push!(Time, solution["time"])
    push!(Arc_reduc, solution["arc_reduc"])
    push!(Node_reduc, solution["node_reduc"])
    push!(Nbr_mtn, solution["nbr_mtn"])
    push!(list_nbr_sts_used, solution["nbr_sts_used"])
    if solution["status"] == MOI.OPTIMAL
        push!(Opt, 1)
    else
        push!(Opt, 0)
    end
    print_solution(solution, Output_file, silent)
     
    dataframe = DataFrames.DataFrame(Instances = Instances, Arc_reduced = Arc_reduc, Node_reduced = Node_reduc, 
                                    Rho = Obj_rho, Lambda = Obj_lambda, Phi = Obj_phi, UB = Obj, LB = Dual_obj, 
                                    Gap = Gap, Nodes = Nbr_nodes, Time = Time, Opt = Opt, 
                                    Nbr_mtn = Nbr_mtn, Nbr_sts_used = list_nbr_sts_used)
    if preprocess
        XLSX.writetable(Outputs_fold*"result_"*inst_name*"_with_prep.xlsx", dataframe, overwrite=true)
    else
        XLSX.writetable(Outputs_fold*"result_"*inst_name*"_without_prep.xlsx", dataframe, overwrite=true)
    end
end 