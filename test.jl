using JuMP
using Gurobi
using LinearAlgebra

"""
Exemple de décomposition de Benders avec callbacks.
Problème : Min c'x + d'y
           s.t. Ax + By >= b
                x >= 0, y >= 0
                x entier

On décompose en:
- Master: Min θ + c'x, s.t. x entier, x >= 0, θ >= lower_bound
- Subproblem: Min d'y, s.t. By >= b - Ax̄, y >= 0
"""

function resoudre_avec_benders()
    # Données du problème
    c = [5.0, 4.0]  # Coûts variables de premier niveau
    d = [6.0, 3.0, 2.0]  # Coûts variables de second niveau
    
    # Contraintes: A*x + B*y >= b
    A = [2.0 1.0;
         1.0 2.0;
         3.0 1.0]
    
    B = [1.0 2.0 1.0;
         2.0 1.0 1.0;
         1.0 1.0 2.0]
    
    b = [10.0, 12.0, 15.0]
    
    # Créer le problème maître
    master = Model(Gurobi.Optimizer)
    set_optimizer_attribute(master, "OutputFlag", 0)
    
    @variable(master, x[1:2] >= 0, Int)
    @variable(master, θ >= -1000)  # Variable pour l'objectif du sous-problème
    @objective(master, Min, dot(c, x) + θ)
    
    # Compteurs
    nb_coupes = Ref(0)
    nb_iterations = Ref(0)
    
    # Fonction pour résoudre le sous-problème
    function resoudre_subproblem(x_val)
        sub = Model(Gurobi.Optimizer)
        set_optimizer_attribute(sub, "OutputFlag", 0)
        
        @variable(sub, y[1:3] >= 0)
        @objective(sub, Min, dot(d, y))
        
        # Contraintes avec x fixé
        rhs = b - A * x_val
        @constraint(sub, con, B * y .>= rhs)
        
        optimize!(sub)
        
        if termination_status(sub) == MOI.OPTIMAL
            return (
                optimal = true,
                obj = objective_value(sub),
                y_val = value.(y),
                dual = dual.(con)
            )
        else
            # Sous-problème infaisable - coupe de faisabilité
            # On récupère le rayon extrême
            return (
                optimal = false,
                dual = dual.(con)
            )
        end
    end
    
    # Callback pour ajouter les coupes de Benders
    function benders_callback(cb_data)
        status = callback_node_status(cb_data, master)
        
        # On n'ajoute des coupes que pour les solutions entières
        if status != MOI.CALLBACK_NODE_STATUS_INTEGER
            return
        end
        
        nb_iterations[] += 1
        
        # Récupérer la solution courante du master
        x_val = callback_value.(cb_data, x)
        θ_val = callback_value(cb_data, θ)
        
        # Résoudre le sous-problème
        result = resoudre_subproblem(x_val)
        
        if result.optimal
            # Coupe d'optimalité
            # θ >= d'y* + π'(b - Ax)
            # où π sont les variables duales
            
            gap = result.obj - θ_val
            
            if gap > 1e-6  # Tolérance
                nb_coupes[] += 1
                π = result.dual
                
                # Construire la coupe: θ >= d'y* + π'(b - Ax)
                # Réarrangé: θ + π'A*x >= π'b
                coupe = @build_constraint(θ + dot(π, A * x) >= dot(π, b))
                
                MOI.submit(master, MOI.LazyConstraint(cb_data), coupe)
                
                println("Iteration $(nb_iterations[]): Coupe #$(nb_coupes[]) ajoutée")
                println("  x = $(round.(x_val, digits=2))")
                println("  θ actuel = $(round(θ_val, digits=2)), obj sous-prob = $(round(result.obj, digits=2))")
                println("  Gap = $(round(gap, digits=2))")
            end
        else
            # Coupe de faisabilité (sous-problème infaisable)
            # π'(b - Ax) <= 0
            nb_coupes[] += 1
            π = result.dual
            
            coupe = @build_constraint(dot(π, A * x) >= dot(π, b))
            MOI.submit(master, MOI.LazyConstraint(cb_data), coupe)
            
            println("Iteration $(nb_iterations[]): Coupe de faisabilité #$(nb_coupes[]) ajoutée")
        end
    end
    
    # Enregistrer le callback
    set_attribute(master, MOI.LazyConstraintCallback(), benders_callback)
    
    # Paramètres Gurobi pour forcer l'utilisation du callback
    set_optimizer_attribute(master, "LazyConstraints", 1)
    set_optimizer_attribute(master, "PreCrush", 1)
    
    println("=== Début de la décomposition de Benders ===\n")
    
    # Résoudre
    optimize!(master)
    
    println("\n=== Solution finale ===")
    println("Statut: $(termination_status(master))")
    println("x = $(value.(x))")
    println("θ = $(value(θ))")
    println("Objectif total = $(objective_value(master))")
    println("Nombre de coupes ajoutées: $(nb_coupes[])")
    println("Nombre d'itérations: $(nb_iterations[])")
    
    # Vérifier avec le sous-problème final
    x_final = value.(x)
    result = resoudre_subproblem(x_final)
    println("\nVérification:")
    println("  Coût premier niveau (c'x) = $(dot(c, x_final))")
    println("  Coût second niveau (d'y) = $(result.obj)")
    println("  Total = $(dot(c, x_final) + result.obj)")
end

# Exécuter
resoudre_avec_benders()