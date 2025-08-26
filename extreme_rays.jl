using JuMP
using Gurobi

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "InfUnbdInfo", 1)

@variable(model, x >= 0)
@constraint(model, c1, x >= 1)
@constraint(model, c2, x <= 0)  # --> modèle infaisable

@objective(model, Min, x)

optimize!(model)

status = termination_status(model)

if status == MOI.INFEASIBLE
    println("Modèle infaisable. Extraction du rayon de Farkas...")

    # Extraire les multiplicateurs de Farkas pour les contraintes
    ray_c1 = shadow_price(c1)
    ray_c2 = shadow_price(c2)

    println("Rayon extrême (Farkas multipliers):")
    println("  λ₁ (pour c1: x ≥ 1) = ", ray_c1)
    println("  λ₂ (pour c2: x ≤ 0) = ", ray_c2)
else
    println("Modèle faisable. Status: ", status)
end
