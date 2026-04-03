using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob


include("model_amrp.jl")
include("../build_graph.jl")
include("functions.jl")
include("rf_and_fo.jl")
env = Gurobi.Env()

#= inst_name = ARGS[1]
preprocess = parse(Bool, ARGS[2])
 =# 
instance = "354FL_8A_"
inst_list = ["136FL_5A_2"]
preprocess = false

nbr_thread = 8
graph_reduc = false
fold_1 = "RESULTS/A_MTN_1/"
Outputs_fold = fold_1

#Pré-run
file_essai = "../INSTANCES/instances_json/A_MTN_2/127FL_5A_1.json"
inst_data = build_graph(file_essai)
path_f = Outputs_fold*"result_127FL_5A_1_RF.txt"
Otp_file = open(path_f, "w")
t_limit = 5
mdl = buildM(env, inst_data, nbr_thread, true, graph_reduc, t_limit)
sz_day, ovlap = 3, 2
sol = relax_and_fix(mdl, inst_data, Otp_file, sz_day, ovlap)

Instances, Obj, Dual_obj, Gap, Nbr_nodes  = [], [], [], [], []
Time, Opt, Nbr_mtn = [], [], []
Nbr_sts_used, list_nbr_sts_used = [], []
RF_obj, RF_time, CM_obj, CM_time = [], [], [], []
for v in 2:2
    inst_name = instance*string(v)
    inst_path = "../INSTANCES/instances_literature_json/A_MTN_1/"*inst_name
    if graph_reduc
        path_file = Outputs_fold * "result_" * inst_name * "_WS_GR.txt"
    else
        path_file = Outputs_fold * "result_" * inst_name * "_WS.txt"
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
    rf_sp_limit = 120
    silent = false
    model = buildM(env, instance_data, nbr_thread, silent, graph_reduc, rf_sp_limit)
    size_day, overlap = 3, 2
    println("RELAX AND FIX")
    solution = relax_and_fix(model, instance_data, Output_file, size_day, overlap)
    rf_obj = solution.obj
    rf_time = solution.other_info["time"]
    push!(RF_obj, rf_obj)
    push!(RF_time, rf_time)
    println("TIME  : ", rf_time, " seconds")

    if rf_obj > 0.0
        write_both(Output_file, "\nCOMPLETE MODEL")
        silent = false
        time_limit = 300
        x_start = solution.x
        y_start = solution.y
        obj_start = solution.obj
        close(Output_file)
        solution = model_amrp_with_warm_start(env, instance_data, x_start, y_start, obj_start, path_file, nbr_thread, silent, graph_reduc, time_limit)
        Output_file = open(path_file, "a") 
        other_info = solution.other_info   
        cm_obj =  solution.obj
        cm_time = other_info["time"]
        push!(CM_obj, cm_obj)
        push!(CM_time, cm_time)
        push!(Dual_obj, other_info["dual_obj"])
        push!(Gap, other_info["gap"])
        push!(Nbr_nodes, other_info["nbr_nodes"])
        push!(Nbr_mtn, other_info["nbr_mtn"])
        push!(Nbr_sts_used, other_info["nbr_sts_used"])

        if other_info["status"] == MOI.OPTIMAL
            push!(Opt, 1)
        else
            push!(Opt, 0)
        end  
    else
        write_both(Output_file, "\nNO COMPLETE MODEL")
        cm_obj = solution.obj
        cm_time = 0.0
        push!(CM_obj, cm_obj)
        push!(CM_time, cm_time)
        push!(Dual_obj, 0.0)
        push!(Gap, 0.0)
        push!(Nbr_nodes, 0.0)
        push!(Nbr_mtn, 0.0)
        push!(Nbr_sts_used, 0.0)
        push!(Opt, 1)
    end 
    write_both(Output_file, "\tCM OBJECTIVE VALUE : $cm_obj")
    write_both(Output_file, "\tCM TIME  : $cm_time seconds")
    model = nothing
end
dataframe = DataFrames.DataFrame(Instances = Instances, RF_Obj = RF_obj, RF_Time = RF_time, CM_Obj = CM_obj, CM_Time = CM_time, 
                                LB = Dual_obj, Gap = Gap, Nodes = Nbr_nodes, Opt = Opt, Nbr_mtn = Nbr_mtn, Nbr_sts_used = Nbr_sts_used)

if graph_reduc
    csv_summary_path = Outputs_fold * "summary_"*instance*"WS_GR.csv"
else
    csv_summary_path = Outputs_fold * "summary_"*instance*"WS.csv"
end
# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

println("\nCSV Summary saved to: ", csv_summary_path)
