using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob

include("__model__ws_milp.jl")
include("../../build_graph.jl")
include("../__functions__.jl")
include("../__rffo__.jl")

env = Gurobi.Env()

nbr_thread = 8
silent = true
graph_reduc = false
all_time_limit = 18000

# =============== Sur Nibi ===============
#instance = ARGS[1]

# =============== En local sur mon ordi ===============
inst = ARGS[1]
instance = "../../INSTANCES/instances_literature_json/A_MTN_1/"*inst

# ============== Gestion des fichiers et la sortie =================
parts = split(instance, "/")
f_fold = parts[end-1]
Output_fold = "RESULTS_WS/" * f_fold * "/"

# ============== Pré-run ================
file_essai = "../../INSTANCES/instances_json/A_MTN_2/127FL_5A_1.json"
inst_data = build_graph(file_essai)
path_f = "RESULTS_WS/result_127FL_5A_1_RF.txt"
Otp_file = open(path_f, "w")
t_limit = 15
mdl = buildM(env, inst_data, nbr_thread, silent, graph_reduc, t_limit)
sz_day, ovlap = 3, 2
sol = relax_and_fix(mdl, inst_data, Otp_file, sz_day, ovlap)

# ============= Run officiel ============== 
inst_name = splitext(basename(instance))[1]

if graph_reduc
    path_file = Output_fold*"result_"*inst_name*"_WS_GR.txt"
else
    path_file = Output_fold*"result_"*inst_name*"_WS.txt"
end 

Output_file = open(path_file, "w")

write_both(Output_file, "======================= INSTANCE "*inst_name*"=======================")
write_both(Output_file, "\nGRAPH RÉDUCTION : "*string(graph_reduc)*"\n")

instance_data = build_graph(instance)
rf_sp_time_limit = 60
model = buildM(env, instance_data, nbr_thread, silent, graph_reduc, rf_sp_time_limit)
size_day, overlap = 3, 2
write_both(Output_file, "RELAX AND FIX")
solution = relax_and_fix(model, instance_data, Output_file, size_day, overlap)
rf_obj = solution.obj
rf_time = solution.other_info["time"]
write_both(Output_file, "\tRF OBJECTIVE VALUE : $rf_obj")
write_both(Output_file, "\tRF TIME  :  $rf_time seconds")

if rf_obj > 0.0
    write_both(Output_file, "\nCOMPLETE MODEL")
    silent = true
    milp_time_limit = all_time_limit-rf_time
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


dataframe = DataFrames.DataFrame(Instances = [inst_name], RF_Obj = [rf_obj], RF_Time = [rf_time], CM_Obj = [cm_obj], CM_Time = [cm_time], 
                                LB = [Dual_obj], Gap = [Gap], Nodes = [Nbr_nodes], Opt = [Opt], Nbr_mtn = [Nbr_mtn], Nbr_sts_used = [Nbr_sts_used])

if graph_reduc
    csv_summary_path = Output_fold * "summary_" * inst_name* "_WS_GR.csv"
else
    csv_summary_path = Output_fold * "summary_" * inst_name* "_WS.csv"
end 
# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

write_both(Output_file, "\nCSV Summary saved to: $csv_summary_path")
