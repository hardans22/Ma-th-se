using JuMP, Gurobi, JSON, CSV, MathOptInterface, DataFrames, XLSX, Glob


include("rf_and_fo.jl")
include("../build_graph.jl")
include("functions.jl")
env = Gurobi.Env()

instance = "354FL_8A_"
inst_list = ["136FL_5A_2"]

nbr_thread = 8
silent = true
time_limit = 15
graph_reduc = false
#model_amrp(file*".json")
fold_1 = "RESULTS/A_MTN_5/"
Outputs_fold = fold_1

#Pré-run
file_essai = "../INSTANCES/instances_json/A_MTN_2/127FL_5A_1.json"
inst_data = build_graph(file_essai)
path_f = Outputs_fold*"result_127FL_5A_1_RF.txt"
Otp_file = open(path_f, "w")
mdl = buildM(env, inst_data, nbr_thread, true, graph_reduc, time_limit)
sz_day, ovlap = 3, 1 
sol = relax_and_fix(mdl, inst_data, Otp_file, sz_day, ovlap)
sol_ = fix_and_optimize(mdl, sol.x, sol.y, inst_data, Otp_file, sz_day, ovlap)

Instances, RF_Obj, FO_Obj  = [], [], []
RF_Time, FO_Time = [], []
list_nbr_sts_used = []

for v in 1:1
    inst_name = instance*string(v)
    inst_path = "../INSTANCES/instances_literature_json/A_MTN_1/"*inst_name
    
    if graph_reduc
        path_file = Outputs_fold * "result_" * inst_name * "_RF_GR.txt"
    else
        path_file = Outputs_fold * "result_" * inst_name * "_RF.txt"
    end
    
    Output_file = open(path_file, "w")
    write_both(Output_file, "======================= INSTANCE "*inst_name*"=======================")
    write_both(Output_file, "\nGRAPH RÉDUCTION : "*string(graph_reduc)*"\n")

    push!(Instances, inst_name)

    instance_file = inst_path*".json"
    instance_data = build_graph(instance_file)
    rf_sp_limit = 60
    model = buildM(env, instance_data, nbr_thread, silent, graph_reduc, rf_sp_limit)
    size_day, overlap = 3, 2
    write_both(Output_file, "RELAX AND FIX")
    solution = relax_and_fix(model, instance_data, Output_file, size_day, overlap)
    rf_obj = solution.obj
    rf_time = solution.other_info["time"]
    push!(RF_Obj, rf_obj)
    push!(RF_Time, rf_time)
    write_both(Output_file, "\tRF OBJECTIVE VALUE : $rf_obj")
    write_both(Output_file, "\tRF TIME  :  $rf_time seconds")

    if solution.obj > 0
        write_both(Output_file, "\nFIX AND OPTIMIZE")
        sx = solution.x 
        sy = solution.y
        model = solution.other_info["model"]
        fo_sp_limit = 60
        #= set_optimizer_attribute(model, "TimeLimit", fo_sp_limit) 
        set_optimizer_attribute(model, "OutputFlag", 1) 
         =#
        size_day, overlap = 3, 2
        solution = fix_and_optimize(model, sx, sy, instance_data, Output_file, size_day, overlap)
        
        fo_obj = solution.obj
        fo_time = solution.other_info["time"]
    else
        write_both(Output_file, "\n NO FIX AND OPTIMIZE")
        fo_obj = solution.obj
        fo_time = 0.0
    end     
    
    push!(FO_Obj, fo_obj)
    push!(FO_Time, fo_time)
    write_both(Output_file, "\tFO OBJECTIVE VALUE : $fo_obj")
    write_both(Output_file, "\tFO TIME  : $fo_time seconds")
    #print_solution(solution, instance_data, Output_file, silent)
    model = nothing
    #close(Output_file)
end 

dataframe = DataFrames.DataFrame(Instances = Instances, RF_Obj = RF_Obj, RF_Time = RF_Time, FO_Obj = FO_Obj, FO_Time = FO_Time)

csv_summary_path = Outputs_fold * "summary_" * instance* "RF.csv"


# --- Exporter les DataFrames en CSV ---
CSV.write(csv_summary_path, dataframe)

println("\nCSV Summary saved to: ", csv_summary_path)
