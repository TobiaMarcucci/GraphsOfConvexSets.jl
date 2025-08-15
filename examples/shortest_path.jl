# Inspired from https://github.com/TobiaMarcucci/gcspy/blob/main/examples/shortest_path/shortest_path.py

using JuMP
import Pajarito, HiGHS, Hypatia
using LinearAlgebra
using GraphsOfConvexSets

model = Model(() -> Optimizer(
    optimizer_with_attributes(
        Pajarito.Optimizer,
        "oa_solver" => optimizer_with_attributes(
            HiGHS.Optimizer,
            MOI.Silent() => false,
            "mip_feasibility_tolerance" => 1e-8,
            "mip_rel_gap" => 1e-6,
        ),
        "conic_solver" =>
            optimizer_with_attributes(Hypatia.Optimizer, MOI.Silent() => false),
    )
))

@variable(model, x[1:5, 1:2])
MOI.set.(model, VariableVertexOrEdge(), x, 1:5)

# Centers
C = [
     1    0
    10    0
     4    2
     5.5 -2
     7    2
]

# source vertex
D = Diagonal([1, 1/2]) # scaling matrix
con_ref = @constraint(model, [1; D * (x[1, :] - C[1, :])] in SecondOrderCone())
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 1)

# target vertex
D = Diagonal([1/2, 1]) # scaling matrix
con_ref = @constraint(model, [1; D * (x[2, :] - C[2, :])] in SecondOrderCone())
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 2)
con_ref = @constraint(model, x[2, 1] <= C[2, 1]) # cut right half of the set
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 2)

# vertex 1
con_ref = @constraint(model, [1; x[3, :] - C[3, :]] in MOI.NormInfinityCone(3))
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 3)

# vertex 2
con_ref = @constraint(model, [1.2; x[4, :] - C[4, :]] in MOI.NormOneCone(3))
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 4)
con_ref = @constraint(model, [1; x[4, :] - C[4, :]] in SecondOrderCone())
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 4)

# vertex 3
con_ref = @constraint(model, [1; x[5, :] - C[5, :]] in SecondOrderCone())
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 5)

edges = [(1, 3), (1, 4), (3, 4), (3, 5), (4, 5), (4, 2), (5, 2)]

for (src, dst) in edges
    # Cost of the edge
    cost = @variable(model)
    MOI.set(model, VariableVertexOrEdge(), cost, (src,dst))
    MOI.set(model, ObjectiveVertexOrEdge((src, dst)), MOI.VariableIndex(cost))

    cons_ref = @constraint(model, [cost; x[dst, :] - x[src, :]] in SecondOrderCone())
    MOI.set(model, ConstraintVertexOrEdge(), cons_ref, (src, dst))

    # Constraint on the edge
    cons_ref = @constraint(model, x[src, 2] <= x[dst, 2])
    MOI.set(model, ConstraintVertexOrEdge(), cons_ref, (src, dst))
end

MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
MOI.set(model, Problem(), ShortestPathProblem(1, 2))

optimize!(model)

sol = MOI.get(model, GraphsOfConvexSets.SubGraph())

using GraphPlot
gplot(sol, nodelabel = 1:5, layout = shell_layout)

xv = value.(x)

using Plots
p = plot()

using Graphs
for edge in edges
    if Graphs.has_edge(sol, edge)
        v = [edge...]
        plot!(p, xv[v, 1], xv[v, 2], arrow = true, label = "", color = :grey, linewidth=2)
        # Is there some way to verify the correctness of the solution (in particular the costs)?
        # Since in formulation 5.8, if the cost is a variable x_v on a vertex, we may not have have x_v* = z_v*
        # Indeed, the set to which x_v belong is {x | x ≥ 0} is unbounded, and his recession cone is {(x,y) | x ≥ 0, y ≥ 0}. Therefore, when y = 0, we may not have x = 0
        # 5.8c gives for y_v = 0, (z_v*,0) ∈ hom(X_v) and (x_v* - z_v*) ∈ X_v. As we want to min z_v*, it will be set to 0, and we have x_v* ∈ X_v
        #            for y_v = 1, z_v* ∈ X_v and (x_v* - z_v*, 0) ∈ X_v. And we may have that x_v* ≠ z_v* as long as x_v* ≥ z_v*

        # If the cost is a variable x_e on an edge, no constraint ensure the equality between x_e* and z_e*
    end
end

scatter!(p, [C[:, 1]], [C[:, 2]], label = "", markerstrokewidth = 0, markersize = 3)

for i in 1:5
    plot!(p, [C[i, 1], xv[i, 1]], [C[i, 2], xv[i, 2]], label = "", linestyle = :dash, color = :lightgrey)
    scatter!(p, [xv[i, 1]], [xv[i, 2]], label = "$i", markerstrokewidth=0)
end

p