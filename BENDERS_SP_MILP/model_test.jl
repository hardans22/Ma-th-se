using JuMP

function model_test(env, instance_data, nbr_thread, silent, time_limit)
    if graph_reduc
        instance_data = graph_reduction(instance_data)
    end

    graph = instance_data.graph
    node_sets = graph.node_sets
    arc_sets = graph.arc_sets
    fl_data = instance_data.fl_data
    other_data = instance_data.other_data

    d = fl_data.d
    H_aircraft = node_sets.H_aircraft
    nbr_K = instance_data.nbr_ac
    fl_nodes = node_sets.fl_nodes
    a_nodes = node_sets.ac_nodes
    L_M = node_sets.mtn_nodes
    V = graph.nodes
    V_wt_st = vcat(a_nodes, fl_nodes) 

    A = graph.arcs
    A_S = arc_sets.arcs_S
    A_T = arc_sets.arcs_T
    A_M = arc_sets.arcs_M
    A_M_bar = arc_sets.arcs_M_bar
    
    if graph_reduc
        F_bar = other_data["F_bar"]
        d_bar = other_data["d_bar"]
    end

    pred_A, succ_A = Dict(), Dict()
    pred_A, succ_A = update_neighborhood(A, pred_A, succ_A)
    pred_A_M, succ_A_M = Dict(), Dict()
    pred_A_M, succ_A_M = update_neighborhood(A_M, pred_A_M, succ_A_M)
    pred_A_M_bar, succ_A_M_bar = Dict(), Dict()
    pred_A_M_bar, succ_A_M_bar = update_neighborhood(A_M_bar, pred_A_M_bar, succ_A_M_bar)
    
    # ===================== Master problem =====================
    
    master_model = Model(() -> Gurobi.Optimizer(env))
    GRBsetintparam(env, "OutputFlag", silent ? 0 : 1)
    set_optimizer_attribute(master_model, "Threads", nbr_thread)
    set_optimizer_attribute(master_model, "TimeLimit", time_limit) 
    set_optimizer_attribute(master_model, "LazyConstraints", 1)
    set_optimizer_attribute(master_model, "IntFeasTol", 1e-9)
    #set_optimizer_attribute(master_model, "Heuristics", 0)

    # Variables
    @variable(master_model, x[k in A], Bin)
    @variable(master_model, Z >= 0)

    # Contraintes
    #= @constraint(master_model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(master_model, sum(x[i] for i in A_T) == nbr_K)
     =#
    @constraint(master_model, mp_c1[i in V_wt_st], sum(x[(i,j)] for j in get(succ_A, i, [])) == 1)
    @constraint(master_model, mp_c2[i in V_wt_st], sum(x[(j,i)] for j in get(pred_A, i, [])) == sum(x[(i,j)] for j in get(succ_A, i, [])))
    
    @constraint(master_model, 210*x[("s", "N121")] + 210*x[("N121", "MLX_ESB_31920_31995")] + 210*x[("MLX_ESB_31920_31995", "ESB_TZX_32040_32120")] + 210*x[("ESB_TZX_32040_32120", "TZX_ESB_32155_32235")] + 210*x[("TZX_ESB_32155_32235", "ESB_SAW_32400_32460")] + 210*x[("ESB_SAW_32400_32460", "SAW_ESB_32520_32580")] + 210*x[("SAW_ESB_32520_32580", "ESB_SZF_32660_32720")] + 210*x[("ESB_SZF_32660_32720", "SZF_ESB_32755_32815")] + 210*x[("SZF_ESB_32755_32815", "ESB_BAL_33505_33595")] - Z <= 1680.0)

    @constraint(master_model, x[("s", "N121")] + x[("N121", "MLX_ESB_31920_31995")] + x[("MLX_ESB_31920_31995", "ESB_SAW_32040_32100")] + x[("ESB_SAW_32040_32100", "SAW_ESB_32160_32220")] + x[("SAW_ESB_32160_32220", "ESB_MQM_32330_32425")] + x[("ESB_MQM_32330_32425", "MQM_ESB_32460_32555")] + x[("MQM_ESB_32460_32555", "ESB_ADA_32640_32700")] + x[("ESB_ADA_32640_32700", "ADA_ESB_32740_32805")] + x[("ADA_ESB_32740_32805", "ESB_MLX_32860_32935")] + x[("ESB_MLX_32860_32935", "MLX_ESB_33360_33435")] + x[("MLX_ESB_33360_33435", "ESB_AJI_33490_33585")] <= 10.0)

    @constraint(master_model, x[("s", "N121")] + x[("N121", "MLX_ESB_31920_31995")] + x[("MLX_ESB_31920_31995", "ESB_TZX_32040_32120")] + x[("ESB_TZX_32040_32120", "TZX_ESB_32155_32235")] + x[("TZX_ESB_32155_32235", "ESB_MQM_32330_32425")] + x[("ESB_MQM_32330_32425", "MQM_ESB_32460_32555")] + x[("MQM_ESB_32460_32555", "ESB_ADA_32640_32700")] + x[("ESB_ADA_32640_32700", "ADA_ESB_32740_32805")] + x[("ADA_ESB_32740_32805", "ESB_MLX_32860_32935")] + x[("ESB_MLX_32860_32935", "MLX_ESB_33360_33435")] <= 9.0)

    @constraint(master_model, x[("s", "N121")] + x[("N121", "MLX_ESB_31920_31995")] + x[("MLX_ESB_31920_31995", "ESB_AJI_32050_32145")] + x[("ESB_AJI_32050_32145", "AJI_ESB_32180_32285")] + x[("AJI_ESB_32180_32285", "ESB_ECN_32340_32410")] + x[("ESB_ECN_32340_32410", "ECN_HTY_32520_32585")] + x[("ECN_HTY_32520_32585", "HTY_ECN_32625_32690")] + x[("HTY_ECN_32625_32690", "ECN_ESB_32735_32805")] + x[("ECN_ESB_32735_32805", "ESB_ECN_32865_32935")] + x[("ESB_ECN_32865_32935", "ECN_ESB_33360_33430")] + x[("ECN_ESB_33360_33430", "ESB_VAN_33565_33660")] <= 10.0)

    @constraint(master_model, x[("s", "N121")] + x[("N121", "MLX_ESB_31920_31995")] + x[("MLX_ESB_31920_31995", "ESB_TZX_32040_32120")] + x[("ESB_TZX_32040_32120", "TZX_ESB_32155_32235")] + x[("TZX_ESB_32155_32235", "ESB_ECN_32340_32410")] + x[("ESB_ECN_32340_32410", "ECN_HTY_32520_32585")] + x[("ECN_HTY_32520_32585", "HTY_ECN_32625_32690")] + x[("HTY_ECN_32625_32690", "ECN_ESB_32735_32805")] + x[("ECN_ESB_32735_32805", "ESB_ECN_32865_32935")] + x[("ESB_ECN_32865_32935", "ECN_ESB_33360_33430")] + x[("ECN_ESB_33360_33430", "ESB_VAN_33565_33660")] <= 10.0)

    @constraint(master_model, x[("s", "N121")] + x[("N121", "MLX_ESB_31920_31995")] + x[("MLX_ESB_31920_31995", "ESB_TZX_32040_32120")] + x[("ESB_TZX_32040_32120", "TZX_ESB_32155_32235")] + x[("TZX_ESB_32155_32235", "ESB_SAW_32400_32460")] + x[("ESB_SAW_32400_32460", "SAW_ESB_32520_32580")] + x[("SAW_ESB_32520_32580", "ESB_SZF_32660_32720")] + x[("ESB_SZF_32660_32720", "SZF_ESB_32755_32815")] + x[("SZF_ESB_32755_32815", "ESB_MLX_32860_32935")] + x[("ESB_MLX_32860_32935", "MLX_ESB_33360_33435")] + x[("MLX_ESB_33360_33435", "ESB_AJI_33490_33585")] <= 10.0)

    @constraint(master_model, x[("s", "N121")] + x[("N121", "MLX_ESB_31920_31995")] + x[("MLX_ESB_31920_31995", "ESB_TZX_32040_32120")] + x[("ESB_TZX_32040_32120", "TZX_ESB_32155_32235")] + x[("TZX_ESB_32155_32235", "ESB_SAW_32400_32460")] + x[("ESB_SAW_32400_32460", "SAW_ESB_32520_32580")] + x[("SAW_ESB_32520_32580", "ESB_SZF_32660_32720")] + x[("ESB_SZF_32660_32720", "SZF_ESB_32755_32815")] + x[("SZF_ESB_32755_32815", "ESB_YEI_32875_32935")] + x[("ESB_YEI_32875_32935", "YEI_ESB_33360_33420")] + x[("YEI_ESB_33360_33420", "ESB_TZX_33480_33560")] + x[("ESB_TZX_33480_33560", "TZX_ESB_33595_33675")] <= 11.0)

    @constraint(master_model, x[("s", "N121")] + x[("N121", "MLX_ESB_31920_31995")] + x[("MLX_ESB_31920_31995", "ESB_TZX_32040_32120")] + x[("ESB_TZX_32040_32120", "TZX_ESB_32155_32235")] + x[("TZX_ESB_32155_32235", "ESB_SAW_32400_32460")] + x[("ESB_SAW_32400_32460", "SAW_ESB_32520_32580")] + x[("SAW_ESB_32520_32580", "ESB_SZF_32660_32720")] + x[("ESB_SZF_32660_32720", "SZF_ESB_32755_32815")] + x[("SZF_ESB_32755_32815", "ESB_MLX_32860_32935")] + x[("ESB_MLX_32860_32935", "MLX_ESB_33360_33435")] + x[("MLX_ESB_33360_33435", "ESB_SAW_33480_33540")] + x[("ESB_SAW_33480_33540", "SAW_ESB_33600_33660")] <= 11.0)
    
    @constraint(master_model, x[("s", "N121")] + x[("N121", "MLX_ESB_31920_31995")] + x[("MLX_ESB_31920_31995", "ESB_TZX_32040_32120")] + x[("ESB_TZX_32040_32120", "TZX_ESB_32155_32235")] + x[("TZX_ESB_32155_32235", "ESB_SAW_32400_32460")] + x[("ESB_SAW_32400_32460", "SAW_ESB_32520_32580")] + x[("SAW_ESB_32520_32580", "ESB_ADA_32640_32700")] + x[("ESB_ADA_32640_32700", "ADA_ESB_32740_32805")] + x[("ADA_ESB_32740_32805", "ESB_MLX_32860_32935")] + x[("ESB_MLX_32860_32935", "MLX_ESB_33360_33435")] + x[("MLX_ESB_33360_33435", "ESB_AJI_33490_33585")] <= 10.0)

    @constraint(master_model, 120*x[("s", "N121")] + 120*x[("N121", "MLX_ESB_31920_31995")] + 120*x[("MLX_ESB_31920_31995", "ESB_TZX_32040_32120")] + 120*x[("ESB_TZX_32040_32120", "TZX_ESB_32155_32235")] + 120*x[("TZX_ESB_32155_32235", "ESB_SAW_32280_32340")] + 120*x[("ESB_SAW_32280_32340", "SAW_AYT_32385_32450")] + 120*x[("SAW_AYT_32385_32450", "AYT_SAW_32490_32555")] + 120*x[("AYT_SAW_32490_32555", "SAW_KYA_32670_32740")] + 120*x[("SAW_KYA_32670_32740", "KYA_SAW_32775_32845")] + 120*x[("KYA_SAW_32775_32845", "SAW_ADB_33360_33420")] - Z <= 1080.0)

    @constraint(master_model, 205*x[("s", "N121")] + 205*x[("N121", "MLX_ESB_31920_31995")] + 205*x[("MLX_ESB_31920_31995", "ESB_TZX_32040_32120")] + 205*x[("ESB_TZX_32040_32120", "TZX_ESB_32155_32235")] + 205*x[("TZX_ESB_32155_32235", "ESB_SAW_32400_32460")] + 205*x[("ESB_SAW_32400_32460", "SAW_ESB_32520_32580")] + 205*x[("SAW_ESB_32520_32580", "ESB_ADA_32640_32700")] + 205*x[("ESB_ADA_32640_32700", "ADA_ESB_32740_32805")] + 205*x[("ADA_ESB_32740_32805", "ESB_BAL_33505_33595")] - Z <= 1640.0)

    @constraint(master_model, x[("s", "N121")] + x[("N121", "MLX_ESB_31920_31995")] + x[("MLX_ESB_31920_31995", "ESB_SAW_32040_32100")] + x[("ESB_SAW_32040_32100", "SAW_ESB_32160_32220")] + x[("SAW_ESB_32160_32220", "ESB_MQM_32330_32425")] + x[("ESB_MQM_32330_32425", "MQM_ESB_32460_32555")] + x[("MQM_ESB_32460_32555", "ESB_SZF_32660_32720")] + x[("ESB_SZF_32660_32720", "SZF_ESB_32755_32815")] + x[("SZF_ESB_32755_32815", "ESB_MLX_32860_32935")] + x[("ESB_MLX_32860_32935", "MLX_ESB_33360_33435")] + x[("MLX_ESB_33360_33435", "ESB_AJI_33490_33585")] <= 10.0)


    @objective(master_model, Min, Z)

    optimize!(master_model)

    println("Z_value = ", value(Z))
end 






