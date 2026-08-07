using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

include("__model__ws_milp.jl")
include("../../build_graph.jl")
include("../__functions__.jl")
include("../BENDERS/__benders__.jl")

env = Gurobi.Env()

nbr_thread = 8
silent = false
graph_reduc = false
all_time_limit = 18000

# =============== Sur Nibi ===============
#instance = ARGS[1]

# =============== En local sur mon ordi ===============
inst = ARGS[1]
instance = "../../INSTANCES/instances_literature_json/"*inst

# ============== Gestion des fichiers et la sortie =================
parts = split(instance, "/")
f_fold = parts[end-1]
Output_fold = "RESULTS_WS/" * f_fold * "/"

# ============= Run officiel ============== 
inst_name = splitext(basename(instance))[1]

if graph_reduc
    path_file = Output_fold*"result_"*inst_name*"_BD_WS_GR.txt"
else
    path_file = Output_fold*"result_"*inst_name*"_BD_WS.txt"
end 

Output_file = open(path_file, "w")

write_both(Output_file, "======================= INSTANCE "*inst_name*"=======================")
write_both(Output_file, "\nGRAPH RÉDUCTION : "*string(graph_reduc)*"\n")

close(Output_file)
sp_method = "HT"
cut_type = "desagg"
instance_data = build_graph(instance)
solution = benders_decomp(env, instance_data, sp_method, cut_type, path_file, nbr_thread, silent, all_time_limit)
other_info = solution.other_info
Output_file = open(path_file, "a") 
#print_solution(solution, instance_data, Output_file)

bd_obj = solution.obj
bd_time = solution.other_info["time"]
write_both(Output_file, "\tBD OBJECTIVE VALUE : $bd_obj")
write_both(Output_file, "\tBD TIME  :  $bd_time seconds")

if bd_obj > 0.0
    write_both(Output_file, "\nCOMPLETE MODEL")
    silent = false
    milp_time_limit = all_time_limit-bd_time
    x_start = solution.x
    y_start = solution.y
    obj_start = solution.obj
    close(Output_file)
    solution = model_amrp_with_warm_start(env, instance_data, x_start, y_start, obj_start, path_file, nbr_thread, silent, graph_reduc, milp_time_limit)
    Output_file = open(path_file, "a") 
    other_info = solution.other_info   
    cm_obj =  solution.obj
    cm_time = other_info["time"]
    Dual_obj = other_info["dual_obj"]
    Gap = other_info["gap"]
    Nbr_nodes = other_info["nbr_nodes"]

    Nbr_mtn = other_info["nbr_mtn"]
    Nbr_sts_used = other_info["nbr_sts_used"]

    if other_info["status"] == MOI.OPTIMAL
        Opt = 1.0
    else
        Opt = 0.0
    end  
else
    write_both(Output_file, "\nNO COMPLETE MODEL")
    cm_obj = solution.obj
    cm_time = 0.0
    Dual_obj = 0.0
    Gap = 0.0
    Nbr_nodes = 0.0 

    Nbr_mtn = 0.0
    Nbr_sts_used = 0.0
    Opt = 1.0
end 
write_both(Output_file, "\t CM OBJECTIVE VALUE : $cm_obj")
write_both(Output_file, "\t CM TIME  : $cm_time seconds")
model = nothing
print_solution(solution, instance_data, Output_file, silent)

#close(Output_file)


dataframe = DataFrames.DataFrame(Instances = [inst_name], BD_Obj = [bd_obj], BD_Time = [bd_time], CM_Obj = [cm_obj], CM_Time = [cm_time], 
                                LB = [Dual_obj], Gap = [Gap], Nodes = [Nbr_nodes], Opt = [Opt], Nbr_mtn = [Nbr_mtn], Nbr_sts_used = [Nbr_sts_used])

if graph_reduc
    csv_summary_path = Output_fold * "summary_" * inst_name* "_BD_WS_GR.csv"
else
    csv_summary_path = Output_fold * "summary_" * inst_name* "_BD_WS.csv"
end 
# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

write_both(Output_file, "\nCSV Summary saved to: $csv_summary_path")
