# Inspired from `gcspy/examples/shortest_path.py`

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

# source vertex
x = @variable(model, [1:2])
MOI.set.(model, VariableVertex(), x, 1)
c = [1, 0] # center of the set
D = Diagonal([1, 1/2]) # scaling matrix
con_ref = @constraint(model, [1; D * (x - c)] in SecondOrderCone())
MOI.set(model, VertexOrEdge(), con_ref, 1)

# target vertex
x = @variable(model, [1:2])
MOI.set.(model, VariableVertex(), x, 2)
c = [10, 0] # center of the set
D = Diagonal([1/2, 1]) # scaling matrix
con_ref = @constraint(model, [1; D * (x - c)] in SecondOrderCone())
MOI.set(model, VertexOrEdge(), con_ref, 2)
con_ref = @constraint(model, x[1] <= c[1]) # cut right half of the set
MOI.set(model, VertexOrEdge(), con_ref, 2)

# vertex 1
x = @variable(model, [1:2])
MOI.set.(model, VariableVertex(), x, 3)
c = [4, 2] # center of the set
con_ref = @constraint(model, [1; x - c] in MOI.NormInfinityCone(3))
MOI.set(model, VertexOrEdge(), con_ref, 3)

# vertex 2
x = @variable(model, [1:2])
MOI.set.(model, VariableVertex(), x, 4)
c = [5.5, -2] # center of the set
con_ref = @constraint(model, [1.2; x - c] in MOI.NormOneCone(3))
MOI.set(model, VertexOrEdge(), con_ref, 4)
con_ref = @constraint(model, [1; x - c] in MOI.NormOneCone(3))
MOI.set(model, VertexOrEdge(), con_ref, 4)

# vertex 3
x = @variable(model, [1:2])
MOI.set.(model, VariableVertex(), x, 5)
c = [7, 2] # center of the set
con_ref = @constraint(model, [1; x - c] in SecondOrderCone())
MOI.set(model, VertexOrEdge(), con_ref, 5)

MOI.set(model, Problem(), ShortestPathProblem(1, 2))

optimize!(model)
