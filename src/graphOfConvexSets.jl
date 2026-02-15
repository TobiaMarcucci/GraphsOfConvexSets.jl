import JuMP
import Graphs

include("MOI_wrapper.jl")

struct GraphOfConvexSets{M <: JuMP.AbstractModel, T, G <: Graphs.AbstractGraph{T}} <: Graphs.AbstractGraph{T}
    model::M
    graph::G

    GraphOfConvexSets(optimizer_factory<:MOI.AbstractOptimizer, g <: Graphs.AbstractGraph) = new(JuMP.Model(GCSOptimizer(optimizer_factory)), g)
end

Graphs.nv(g::GraphOfConvexSets) = Graphs.nv(g.graph)
Graphs.vertices(g::GraphOfConvexSets) = Graphs.vertices(g.graph)
Graphs.has_vertex(g::GraphOfConvexSets, v::Integer) = Graphs.has_vertex(g.graph, v)
Graphs.add_vertex!(g::GraphOfConvexSets) = Graphs.add_vertex!(g.graph)
Graphs.rem_vertex!(g::GraphOfConvexSets, v::Integer) = Graphs.rem_vertex!(g.graph, v)
# TODO : Remove variables, constraints and cost associated to v when removing v

Graphs.ne(g::GraphOfConvexSets) = Graphs.ne(g.graph)
Graphs.edgetype(g::GraphOfConvexSets) = Graphs.edgetype(g.graph)
Graphs.edges(g::GraphOfConvexSets) = Graphs.edges(g.graph)
Graphs.has_edge(g::GraphOfConvexSets, s<:Integer, d<:Integer) = Graphs.has_edge(g.graph, s, d)
Graphs.add_edge!(g::GraphOfConvexSets, s<:Integer, d<:Integer) = Graphs.add_edge!(g.graph, s, d)
Graphs.rem_edge!(g::GraphOfConvexSets, s<:Integer, d<:Integer) = Graphs.rem_edge!(g.graph, s, d)
# TODO : Remove variables, constraints and cost associated to (s,d) when removing (s,d)

Graphs.is_directed(g::GraphOfConvexSets) = Graphs.is_directed(g.graph)
Graphs.inneighbors(g::GraphOfConvexSets, v<:Integer) = Graphs.inneighbors(g.graph, v)
Graphs.outneighbors(g::GraphOfConvexSets, v<:Integer) = Graphs.outneighbors(g.graph, v)

struct Vertex{M <: JuMP.AbstractModel, T<:Integer} <: JuMP.AbstractModel
    model::M
    vertex::T

    function Vertex(g::GraphOfConvexSets{M,T,G}, v::T) where {M <: JuMP.AbstractModel, T<:Integer, G<:Graphs.AbstractGraph{T}}
        Graphs.has_vertex(g,v) || throw(BoundsError("The graph g does not have a vertex $(v). You can add it through the add_vertex! method."))
        return new(g.model, v)
    end
end

function JuMP.add_variable(v::Vertex, var...)
    var_ref = JuMP.add_variable(v.model, var...)
    MOI.set(v.model, VariableVertexOrEdge(), var_ref, v.vertex)
end

function JuMP.add_constraint(v::Vertex, con)
    _check(v.model, con, con, v.vertex)
    con_ref = JuMP.add_constraint(v.model, con)
    MOI.set(v.model, ConstraintVertexOrEdge(), con_ref, v.vertex)
end

function JuMP.set_objective_function(v::Vertex, func)
    _check(v.model, func, func, v.vertex)
    MOI.set(v.model, VertexOrEdgeObjective(v), moi_function(func))
end

struct Edge{M <: JuMP.AbstractModel, T<:Integer} <: JuMP.AbstractModel
    model::M
    edge::Tuple{T, T}

    function Edge(g::GraphOfConvexSets{M,T,G}, s::T, d::T) where {M <: JuMP.AbstractModel, T <: Integer, G<: Graphs.AbstractGraph{T}}
        Graphs.has_edge(g,s,d) || throw(BoundsError("The graph g does not have an edge from $(s) to $(d). You can add it through the add_edge! method."))
    end
end

function JuMP.add_variable(e::Edge, var)
    var_ref = JuMP.add_variable(e.model, var)
    MOI.set(e.model, VariableVertexOrEdge(), var_ref, e.edge)
end

function JuMP.add_constraint(e::Edge, con)
    _check(e.model, con, con, e.edge)
    con_ref = JuMP.add_constraint(e.model, con)
    MOI.set(e.model, ConstraintVertexOrEdge(), con_ref, e.edge)
end

function JuMP.set_objective_function(e::Edge, func)
    _check(e.model, func, func, e.edge)
    MOI.set(e.model, VertexOrEdgeObjective(e), moi_function(func))
end

JuMP.set_objective_sense(g::GraphOfConvexSets, sense::MOI.OptimizationSense) = MOI.set_objective_sense(g.model, sense) # Does it make sense as it is intrisincly relative to the problem, i.e., a shortest path problem will be a min problem?

# TODO : extract model associated to a given vertex/edge (see what has been done for the serializer)

function Graphs.shortest_path(g::GraphOfConvexSets{M,T,G}, s::T, t::T) where {T}
    Graphs.has_vertex(g, s) || Graphs.has_vertex(g, t) || throw(ArgumentError("The graph g does not have either vertex $(s), either vertex $(t)."))

    MOI.set_objective_sense(g.model, sense)
    MOI.set(g.model, Problem(), ShortestPathProblem(1,2))
    optimize!(g.model)
    # TODO
end