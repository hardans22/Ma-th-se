function model_test(env, instance_data, benders_cuts, nbr_thread, silent, time_limit)
    graph       = instance_data.graph
    node_sets   = graph.node_sets
    arc_sets    = graph.arc_sets
    fl_data     = instance_data.fl_data
    other_data  = instance_data.other_data

    

    H_aircraft  = node_sets.H_aircraft
    fl_nodes    = node_sets.fl_nodes
    a_nodes     = node_sets.ac_nodes
    V_wt_st     = vcat(a_nodes, fl_nodes) 

    A   = graph.arcs
    A_S = arc_sets.arcs_S
    A_T = arc_sets.arcs_T
    A_M = arc_sets.arcs_M
    L_M = node_sets.mtn_nodes
    V   = graph.nodes
    d   = fl_data.d

    nbr_K   = instance_data.nbr_ac
    A_M_bar = arc_sets.arcs_M_bar
    
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
    set_optimizer_attribute(master_model, "FeasibilityTol", 1e-09)
    set_optimizer_attribute(master_model, "Cuts", 0)
    # Variables
    master_model[:x] = @variable(master_model, x[arc in A], Bin)
    master_model[:Z] = @variable(master_model, Z[k in a_nodes] >= 0)

    # Contraintes
    #= @constraint(master_model, sum(x[i] for i in A_S) == nbr_K)
    @constraint(master_model, sum(x[i] for i in A_T) == nbr_K)
     =#
    @constraint(master_model, mp_c1[i in V_wt_st], sum(x[(i,j)] for j in get(succ_A, i, [])) == 1)
    @constraint(master_model, mp_c2[i in V_wt_st], sum(x[(j,i)] for j in get(pred_A, i, [])) == sum(x[(i,j)] for j in get(succ_A, i, [])))
    
    @objective(master_model, Min, sum(Z[k] for k in a_nodes))

    parse_and_add_constraints!(master_model, benders_cuts)

    optimize!(master_model)
    
    Z_val = value(Z)
    println(Z_val)
end 


function parse_affine_expr(model, expr_str, x, Z)
    expr = AffExpr(0.0)
    
    # Regex pour capturer chaque terme: [coeff] x[("a","b")] ou [coeff] Z[key]
    term_pattern = r"([+-]?\s*[\d.]*)\s*(x\[\(\"([^\"]+)\",\s*\"([^\"]+)\"\)\]|Z\[(\w+)\])"
    
    for m in eachmatch(term_pattern, expr_str)
        coeff_str = replace(strip(m[1]), " " => "")
        coeff = if coeff_str in ("", "+") 
            1.0
        elseif coeff_str == "-"
            -1.0
        else
            parse(Float64, coeff_str)
        end
        
        if m[5] !== nothing
            # Variable Z[key]
            key = m[5]
            add_to_expression!(expr, coeff, Z[key])
        else
            # Variable x[("from","to")]
            from_node = m[3]
            to_node   = m[4]
            add_to_expression!(expr, coeff, x[(from_node, to_node)])
        end
    end
    
    return expr
end

function parse_and_add_constraints!(model, filepath)
    # Récupérer les variables existantes du modèle
    x = model[:x]
    Z = model[:Z]
    
    open(filepath, "r") do file
        for line in eachline(file)
            line = strip(line)
            isempty(line) && continue
            
            # Extraire la partie entre les premières parenthèses
            # Format: ScalarConstraint{...}(EXPR, LessThan{Float64}(RHS))
            m = match(r"ScalarConstraint\{.*?\}\((.*),\s*MathOptInterface\.LessThan\{Float64\}\(([-\d.]+)\)\)", line)
            m === nothing && continue
            
            expr_str = m[1]
            rhs = parse(Float64, m[2])
            
            # Parser l'expression affine
            expr = parse_affine_expr(model, expr_str, x, Z)
            
            @constraint(model, expr <= rhs)
        end
    end
end
