using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

include("../build_graph.jl")
include("__benders__.jl")
include("__functions__.jl")
include("model_test.jl")
env = Gurobi.Env()
#= inst_name = ARGS[1]
preprocess = parse(Bool, ARGS[2])
 =# 
instance = "667FL_20A_"

nbr_thread = 8
silent = false
time_limit = 20
graph_reduc = false
#model_amrp(file*".json")

s_fold = "A_MTN_8/"
Outputs_fold = "RESULTS_BENDERS/"*s_fold
sp_method = "PD"

Instances, Obj, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], []
Time, Opt, nbr_mtn = [], [], []
nbr_feas_cuts, nbr_opti_cuts, nbr_cuts, nbr_iter = [], [], [], []
MP_time, SP_time = [], []

#for v in [3,4,7,8,10]
for v in 5:5
    inst_name = instance*string(v)
    file_path = "../INSTANCES/instances_literature_json/"*s_fold*inst_name
    
    path_file = Outputs_fold*"result_"*inst_name*"_benders.txt"
    Output_file = open(path_file, "w") 
    write_both(Output_file, "-----------------------INSTANCE $inst_name--------------------------")
    push!(Instances, inst_name)
    close(Output_file)
    instance_file = file_path*".json"
    instance_data = build_graph(instance_file)

    #model_test(env, instance_data, nbr_thread, silent, time_limit)
    cuts_path = "/home/dansou/Documents/PhD_studies/Recherche/Code/BENDERS_METHOD/Benders_cuts.txt"
    model_test(env, instance_data, cuts_path, nbr_thread, silent, time_limit)
end 
