# Inspired from https://github.com/TobiaMarcucci/gcspy/blob/main/examples/shortest_path/shortest_path.py

using JuMP
import Pajarito, HiGHS, Hypatia
using LinearAlgebra
using GraphsOfConvexSets
import Graphs

g = GraphOfConvexSets(optimizer_with_attributes(
        Pajarito.Optimizer,
        "oa_solver" => optimizer_with_attributes(
            HiGHS.Optimizer,
            MOI.Silent() => false,
            "mip_feasibility_tolerance" => 1e-8,
            "mip_rel_gap" => 1e-6,
        ),
        "conic_solver" =>
            optimizer_with_attributes(Hypatia.Optimizer, MOI.Silent() => false),
    ),
    Graphs.SimpleDiGraph()
)

# Construct the graph
Graphs.add_vertices!(g, 5)
edges = [(1, 3), (1, 4), (3, 4), (3, 5), (4, 5), (4, 2), (5, 2)]
for (u,v) in edges Graphs.add_edge!(g, u, v) end

# Create programs on vertices
x = Vector{Vector{VariableRef}}(undef, 0)
for v in Graphs.vertices(g) push!(x, @variable(Vertex(g, v), [1:2])) end

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
@edit @constraint(Vertex(g, 1), D * (x[1][:] - C[1, :]) in SecondOrderCone())

# target vertex
D = Diagonal([1/2, 1]) # scaling matrix
@edit @constraint(Vertex(g, 2), [1; D * (x[2][:] - C[2, :])] in SecondOrderCone())
@edit @constraint(Vertex(g, 2), x[2][1] <= C[2, 1]) # cut right half of the set

# vertex 1
@edit @constraint(Vertex(g, 3), [1; x[3][:] - C[3, :]] in MOI.NormInfinityCone(3))

# vertex 2
@edit @constraint(Vertex(g, 4), [1.2; x[4][:] - C[4, :]] in MOI.NormOneCone(3))
@edit @constraint(Vertex(g, 4), [1; x[4][:] - C[4, :]] in SecondOrderCone())

# vertex 3
@edit @constraint(Vertex(g, 5), [1; x[5][:] - C[5, :]] in SecondOrderCone())

edges = [(1, 3), (1, 4), (3, 4), (3, 5), (4, 5), (4, 2), (5, 2)]

cost = Vector{VariableRef}
for (src, dst) in edges
    edge = Edge(g, src, dst)
    # Cost of the edge
    @edit push!(cost, @variable(edge))
    JuMP.set_objective_function(edge, cost) # TODO : find better
    @constraint(edge, [cost; x[dst, :] - x[src, :]] in SecondOrderCone())

    # Constraint on the edge
    @constraint(edge, x[src, 2] <= x[dst, 2])
end

shortest_path(g, 1, 2)