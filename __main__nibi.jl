using JuMP, Gurobi, JSON, MathOptInterface, DataFrames, XLSX, Glob

#include("model_amrp.jl")
include("model_amrp_new.jl")

include("functions.jl")

instance = ARGS[1]
preprocess = parse(Bool, ARGS[2])
#= instance = "./instances_json/104FL_10A_3D_5MS_1.json"
preprocess = false
 =#
nbr_thread = 8
silent = true
time_limit = 7200
FH = true
FC = false
DY = false
MTN_CAP = false

parts = split(instance, "/")
f_fold = parts[end-1]
Outputs_fold = "RESULTS/"*f_fold*"/"

Instances, Obj, Obj_rho, Obj_lambda, Obj_phi, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], [], [], [], []
Time, Opt, Arc_reduc, Node_reduc, Nbr_mtn = [], [], [], [], []
list_nbr_sts_used, Feasibilities = []

inst_name = splitext(basename(instance))[1]
println(inst_name)
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
write(Output_file, "\n-----------------------INSTANCE "*string(inst_name)*"--------------------------")
push!(Instances, inst_name)
solution = model_amrp(instance, nbr_thread, silent, preprocess, time_limit, FH, FC, DY, MTN_CAP)
println("AVANT CALCUL DES AUTRES INDICATEURS")
println("obj_rho = ", solution["obj_rho"])
println("obj_lambda = ", solution["obj_lambda"])
println("obj_phi = ", solution["obj_phi"])

write(Output_file, "\nAVANT CALCUL DES AUTRES INDICATEURS")
write(Output_file, "\nObj_rho = "*string(solution["obj_rho"]))
write(Output_file, "\nObj_lambda = "*string(solution["obj_lambda"]))
write(Output_file, "\nObj_phi = "*string(solution["obj_phi"]))


solution = compute_indicators(solution, FH, FC, DY)
println("APRÈS CALCUL DES AUTRES INDICATEURS")
println("obj_rho = ", solution["obj_rho"])
println("obj_lambda = ", solution["obj_lambda"])
println("obj_phi = ", solution["obj_phi"])

write(Output_file, "\nAPRÈS CALCUL DES AUTRES INDICATEURS")
write(Output_file, "\nObj_rho = "*string(solution["obj_rho"]))
write(Output_file, "\nObj_lambda = "*string(solution["obj_lambda"]))
write(Output_file, "\nObj_phi = "*string(solution["obj_phi"]))

#solution["feasible"] = 1
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
push!(Feasibilities, solution["feasible"])
if solution["status"] == MOI.OPTIMAL
    push!(Opt, 1)
else
    push!(Opt, 0)
end
print_solution(solution, Output_file, silent)

dataframe = DataFrames.DataFrame(Instances = Instances, Arc_reduced = Arc_reduc, Node_reduced = Node_reduc, 
                            Rho = Obj_rho, Lambda = Obj_lambda, Phi = Obj_phi, UB = Obj, LB = Dual_obj, 
                            Gap = Gap, Nodes = Nbr_nodes, Time = Time, Opt = Opt, 
                            Nbr_mtn = Nbr_mtn, Nbr_sts_used = list_nbr_sts_used, Feasible = Feasibilities)
if preprocess
    XLSX.writetable(Outputs_fold*"result_"*inst_name*"_with_prep.xlsx", dataframe, overwrite=true)
else
    XLSX.writetable(Outputs_fold*"result_"*inst_name*"_without_prep.xlsx", dataframe, overwrite=true)
end