using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob


include("../build_graph.jl")
include("__benders__.jl")
include("functions_for_benders.jl")
include("functions.jl")
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

Outputs_fold = "RESULTS_BENDERS/"

Instances, Obj, Obj_rho, Obj_lambda, Obj_phi, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], [], [], [], []
Time, Opt, Arc_reduc, Node_reduc, Nbr_mtn, Feasibilities = [], [], [], [], [], []
list_nbr_sts_used, prep1_time, prep2_time = [], [], []
nbr_feasibility_cuts, nbr_optimality_cuts, nbr_cuts, nbr_iterations = [], [], [], []
nbr_cover_cuts, MP_time, SP_time = [], [], []
#for v in [3,4,7,8,10]
for v in 1:1
    inst_name = instance*string(v)
    file_path = "../INSTANCES/instances_literature_json/A_MTN_1/"*inst_name
    
    path_file = Outputs_fold*"result_"*inst_name*"_benders.txt"
    Output_file = open(path_file, "w") 
    write_both(Output_file, "-----------------------INSTANCE $inst_name--------------------------")
    push!(Instances, inst_name)
    close(Output_file)
    instance_file = file_path*".json"
    instance_data = build_graph(instance_file)
    
    solution = benders_decomp(env, instance_data, path_file, nbr_thread, silent, time_limit)
    Output_file = open(path_file, "a") 
    
    push!(Obj, solution["obj"])
    push!(MP_time, solution["mp_time"])
    push!(SP_time, solution["sp_time"])
    push!(Time, solution["time"])
    push!(nbr_feasibility_cuts, solution["nbr_feasibility_cuts"])
    push!(nbr_cover_cuts, solution["nbr_cover_cuts"])
    push!(nbr_optimality_cuts, solution["nbr_optimality_cuts"])
    push!(nbr_cuts, solution["nbr_cuts"])
    push!(nbr_iterations, solution["nbr_iterations"])
    write_both(Output_file, "\n=== Résultat final ===")
    write_both(Output_file, "Objectif: $(solution["obj"])")
    write_both(Output_file, "Temps total du master: $(round(solution["mp_time"], digits=2))s")
    write_both(Output_file, "Temps total du sp: $(round(solution["sp_time"], digits=2))s") 
    write_both(Output_file, "Temps total: $(solution["time"])s")
    write_both(Output_file, "Nombre de coupes: $(solution["nbr_cuts"])")
    write_both(Output_file, "Nombre de coupes faisabilité: $(solution["nbr_feasibility_cuts"])")
    write_both(Output_file, "Nombre de cover inequalities: $(solution["nbr_cover_cuts"])")
    write_both(Output_file, "Nombre de coupes d'optimalité: $(solution["nbr_optimality_cuts"])")
    write_both(Output_file, "Nombre d'itérations: $(solution["nbr_iterations"])")
end 
dataframe = DataFrames.DataFrame(Instances = Instances, UB = Obj, MP_Time = MP_time, SP_Time = SP_time, 
                                Time = Time, Nbr_iterations = nbr_iterations, Nbr_feasibility_cuts = nbr_feasibility_cuts,
                                Nbr_optimality_cuts = nbr_optimality_cuts, Nbr_cover_cuts = nbr_cover_cuts, 
                                Nbr_cuts = nbr_cuts)

csv_summary_path = Outputs_fold*"summary_"*instance*"benders.csv"
# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

println("CSV Summary saved to: ", csv_summary_path)


