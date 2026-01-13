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

# source
@variable(model, s[1:2])
MOI.set.(model, VariableVertexOrEdge(), s, 1)
cons_ref = @constraint(model, s .== [-2, -2])
MOI.set.(model, ConstraintVertexOrEdge(), cons_ref, 1)

# affine subspaces
@variable(model, q_0[1:3, 1:2])
MOI.set.(model, VariableVertexOrEdge(), q_0, 2:4)
@variable(model, q_1[1:3, 1:2])
MOI.set.(model, VariableVertexOrEdge(), q_1, 2:4)

# vertex 1 : y = x + 1 and x = -17
con_ref = @constraint(model, q_0[1,2] == q_0[1,1] + 1.)
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 2)
con_ref = @constraint(model, q_0[1,1] == -17.)
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 2)
con_ref = @constraint(model, q_1[1,2] == q_1[1,1] + 1.)
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 2)
con_ref = @constraint(model, q_1[1,1] == -17.)
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 2)

cost = @variable(model)
MOI.set(model, VariableVertexOrEdge(), cost, 2)
set_vertex_or_edge_objective(model, 2, cost)
con_ref = @constraint(model, [cost; 0.5; q_1[1,:] - q_0[1,:]] in RotatedSecondOrderCone())
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 2)

# vertex 2 : y = 0.75x
con_ref = @constraint(model, 0.75q_0[2,1] - q_0[2,2] == 0.)
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 3)
con_ref = @constraint(model, 0.75q_1[2,1] - q_1[2,2] == 0.)
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 3)

cost = @variable(model)
MOI.set(model, VariableVertexOrEdge(), cost, 3)
set_vertex_or_edge_objective(model, 3, cost)
con_ref = @constraint(model, [cost; 0.5; q_1[2,:] - q_0[2,:]] in RotatedSecondOrderCone())
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 3)

# vertex 3 : (x,y) ∈ R^2
# con_ref = @constraint(model, 3q_0[3,1] + q_0[3,2] == 3.)
# MOI.set(model, ConstraintVertexOrEdge(), con_ref, 4)
# con_ref = @constraint(model, 3q_1[3,1] + q_1[3,2] == 3.)
# MOI.set(model, ConstraintVertexOrEdge(), con_ref, 4)

cost = @variable(model)
MOI.set(model, VariableVertexOrEdge(), cost, 4)
set_vertex_or_edge_objective(model, 4, 10.0cost)
con_ref = @constraint(model, [cost; 0.5; q_1[3,:] - q_0[3,:]] in RotatedSecondOrderCone())
MOI.set(model, ConstraintVertexOrEdge(), con_ref, 4)

# target
@variable(model, t[1:2])
MOI.set.(model, VariableVertexOrEdge(), t, 5)
cons_ref = @constraint(model, t .== [0., 0.])
MOI.set.(model, ConstraintVertexOrEdge(), cons_ref, 5)

# Edges with source and target
cons_ref = @constraint(model, s .== q_0[3,:])
MOI.set(model, ConstraintVertexOrEdge(), cons_ref[1], (1, 4))
MOI.set(model, ConstraintVertexOrEdge(), cons_ref[2], (1, 4))
cons_ref = @constraint(model, t .== q_1[2,:])
MOI.set(model, ConstraintVertexOrEdge(), cons_ref[1], (3, 5))
MOI.set(model, ConstraintVertexOrEdge(), cons_ref[2], (3, 5))
cons_ref = @constraint(model, t .== q_1[3,:])
MOI.set(model, ConstraintVertexOrEdge(), cons_ref[1], (4, 5))
MOI.set(model, ConstraintVertexOrEdge(), cons_ref[2], (4, 5))

# Constraints between affine subspaces
edges = [(2, 3), (2,4), (3,2), (3,4), (4,2), (4,3)]

for (src, dst) in edges
    cons_ref = @constraint(model, q_1[src-1,:] .== q_0[dst-1,:])
    MOI.set(model, ConstraintVertexOrEdge(), cons_ref[1], (src, dst))
    MOI.set(model, ConstraintVertexOrEdge(), cons_ref[2], (src, dst))
end

MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
MOI.set(model, Problem(), ShortestPathProblem(1, 5))

optimize!(model)

println(value.(q_0))
println(value.(q_1))