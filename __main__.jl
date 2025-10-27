using JuMP
using Gurobi
using JSON
using MathOptInterface
using DataFrames
using XLSX

include("model_amrp_new.jl")
include("benders_decomp.jl")
include("functions.jl")
#= inst_name = ARGS[1]
preprocess = parse(Bool, ARGS[2])
 =# 
instance = "196FL_7A_"
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
time_limit = 20
FH = true
FC = false
DY = false
MTN_CAP = false
#model_amrp(file*".json")

Outputs_fold = "RESULTS/"


for benders in [false]
    Instances, Obj, Obj_rho, Obj_lambda, Obj_phi, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], [], [], [], []
    Time, Opt, Arc_reduc, Node_reduc, Nbr_mtn, Feasibilities = [], [], [], [], [], []
    list_nbr_sts_used = []
    for v in [1, 8, 10]
        inst_name = instance*string(v)
        file_path = "INSTANCES/instances_json/A_MTN_2/"*inst_name
        if benders
            path_file = Outputs_fold*"result_"*inst_name*"_benders.txt"
        else 
            path_file = Outputs_fold*"result_"*inst_name*".txt"
        end
        Output_file = open(path_file, "w") 
        write(Output_file, "INSTANCE "*inst_name)
        println("\nPREPROCESSING : ", preprocess)
        write(Output_file, "\nPREPROCESSING : "*string(preprocess))
        
        println("-----------------------INSTANCE $inst_name--------------------------")
        #write(Output_file, "\n-----------------------INSTANCE "*inst_name*"--------------------------")
        push!(Instances, inst_name)
        if benders
            solution = benders_decomp1(file_path*".json", nbr_thread, silent, time_limit)
        else
            solution = model_amrp(file_path*".json", nbr_thread, silent, preprocess, time_limit, FH, FC, DY, MTN_CAP)
        end
        #= println("obj_rho = ", solution["obj_rho"])
        println("obj_lambda = ", solution["obj_lambda"])
        println("obj_phi = ", solution["obj_phi"])
         =#
        solution = compute_indicators(solution, FH, FC, DY)
        #= println("obj_rho = ", solution["obj_rho"])
        println("obj_lambda = ", solution["obj_lambda"])
        println("obj_phi = ", solution["obj_phi"])
         =#
        #solution["feasible"] = 1
        push!(Obj, solution["obj"])
        push!(Time, solution["time"])
        push!(Obj_rho, solution["obj_rho"])
        push!(Obj_lambda, solution["obj_lambda"])
        push!(Obj_phi, solution["obj_phi"])
        push!(Dual_obj, solution["dual_obj"])
        push!(Gap, solution["gap"])
        push!(Nbr_nodes, solution["nbr_nodes"])
        push!(Arc_reduc, solution["arc_reduc"])
        push!(Node_reduc, solution["node_reduc"])
        push!(Nbr_mtn, solution["nbr_mtn"])
        push!(list_nbr_sts_used, solution["nbr_sts_used"])
        push!(Feasibilities, solution["feasible"])
        
        if solution["status"] == MOI.OPTIMAL
            push!(Opt, 1)
        else
            push!(Opt, 0)
        end
        
        print_solution(solution, Output_file, silent)
    end 
    dataframe = DataFrames.DataFrame(Instances = Instances, Arc_reduced = Arc_reduc, Node_reduced = Node_reduc, 
                                        Rho = Obj_rho, Lambda = Obj_lambda, Phi = Obj_phi, UB = Obj, LB = Dual_obj, 
                                        Gap = Gap, Nodes = Nbr_nodes, Time = Time, Feasible = Feasibilities, Opt = Opt, 
                                        Nbr_mtn = Nbr_mtn, Nbr_sts_used = list_nbr_sts_used)
    
    #dataframe = DataFrames.DataFrame(Instances = Instances, UB = Obj, Time = Time)
    if benders
        XLSX.writetable(Outputs_fold*"RESULTS_"*instance*"benders.xlsx", dataframe, overwrite=true)
    else
        XLSX.writetable(Outputs_fold*"RESULTS_"*instance*".xlsx", dataframe, overwrite=true)
    end 
end