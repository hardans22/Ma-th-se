using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

include("../build_graph.jl")
include("__benders__.jl")
include("__functions__.jl")
env = Gurobi.Env()
#= inst_name = ARGS[1]
preprocess = parse(Bool, ARGS[2])
 =# 
instance = "354FL_8A_"

nbr_thread = 8
silent = false
time_limit = 30
graph_reduc = false
#model_amrp(file*".json")

s_fold = "A_MTN_5/"
Outputs_fold = "RESULTS_BENDERS/"*s_fold
sp_method = "PD"

Instances, Obj, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], []
Time, Opt, nbr_mtn = [], [], []
nbr_feas_cuts, nbr_opti_cuts, nbr_cuts, nbr_iter = [], [], [], []
MP_time, SP_time = [], []

#for v in [3,4,7,8,10]
for v in 30:30
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
    
    solution = benders_decomp(env, instance_data, sp_method, path_file, nbr_thread, silent, time_limit)
    Output_file = open(path_file, "a") 
    other_info = solution.other_info
    #print_solution(solution, instance_data, Output_file, silent)
    
    push!(Obj, solution.obj)
    push!(MP_time, other_info["mp_time"])
    push!(SP_time, other_info["sp_time"])
    push!(Time, other_info["time"])
    push!(nbr_mtn, other_info["nbr_mtn"])
    push!(nbr_feas_cuts, other_info["nbr_feas_cuts"])
    push!(nbr_opti_cuts, other_info["nbr_opti_cuts"])
    push!(nbr_cuts, other_info["nbr_cuts"])
    push!(nbr_iter, other_info["nbr_iter"])
    write_both(Output_file, "\n=== Résultat final ===")
    write_both(Output_file, "Objectif: $(solution.obj)")
    write_both(Output_file, "Temps total du master: $(round(other_info["mp_time"], digits=4))s")
    write_both(Output_file, "Temps total du sp: $(round(other_info["sp_time"], digits=4))s") 
    write_both(Output_file, "Temps total: $(other_info["time"])s")
    write_both(Output_file, "Nombre de coupes: $(other_info["nbr_cuts"])")
    write_both(Output_file, "Nombre de coupes faisabilité: $(other_info["nbr_feas_cuts"])")
    write_both(Output_file, "Nombre de coupes d'optimalité: $(other_info["nbr_opti_cuts"])")
    write_both(Output_file, "Nombre d'itérations: $(other_info["nbr_iter"])")
end 
dataframe = DataFrames.DataFrame(Instances = Instances, UB = Obj, MP_Time = MP_time, SP_Time = SP_time, 
                                Time = Time, Nbr_mtn = nbr_mtn, Nbr_iter = nbr_iter, Nbr_feasi_cuts = nbr_feas_cuts,
                                Nbr_opti_cuts = nbr_opti_cuts, Nbr_cuts = nbr_cuts)

csv_summary_path = Outputs_fold*"summary_"*instance*"benders_"*sp_method*".csv"
# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

println("CSV Summary saved to: ", csv_summary_path)

