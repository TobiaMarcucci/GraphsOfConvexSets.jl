import Graphs
import MathOptInterface as MOI
using LinearAlgebra

const VertexOrEdgeType = Union{Int,Graphs.Edge{Int}}

struct Problem <: MOI.AbstractModelAttribute end
struct ObjectiveVertexOrEdge <: MOI.AbstractModelAttribute
    vertex_or_edge::Union{Int, Tuple{Int, Int}}
end
struct VariableVertexOrEdge <: MOI.AbstractVariableAttribute end
struct ConstraintVertexOrEdge <: MOI.AbstractConstraintAttribute end

mutable struct Optimizer{O} <: MOI.AbstractOptimizer
    inner::O
    graph::Union{Nothing,Graphs.SimpleDiGraph{Int}}
    variable_vertex_or_edge::Dict{MOI.VariableIndex,VertexOrEdgeType} # x -> vertex or edge
    x::Union{Nothing,Vector{MOI.VariableIndex}}
    y::Dict{VertexOrEdgeType,MOI.VariableIndex}
    z::Dict{Tuple{MOI.VariableIndex,MOI.VariableIndex},MOI.VariableIndex} # (x, y) -> z
    function Optimizer(optimizer_constructor)
        inner = MOI.instantiate(optimizer_constructor, with_bridge_type=Float64)
        return new{typeof(inner)}(
            inner,
            nothing,
            Dict{MOI.VariableIndex,Int}(),
            nothing,
            Dict{VertexOrEdgeType,MOI.VariableIndex}(),
            Dict{Tuple{MOI.VariableIndex,MOI.VariableIndex},MOI.VariableIndex}(),
        )
    end
end

MOI.is_empty(model::Optimizer) = MOI.is_empty(model.inner)
function MOI.empty!(model::Optimizer)
    MOI.empty!(model.inner)
    model.graph = nothing
    empty!(model.variable_vertex_or_edge)
    model.x = nothing
    empty!(model.y)
    empty!(model.z)
    return
end

MOI.optimize!(model::Optimizer) = MOI.optimize!(model.inner)

MOI.supports(model::Optimizer, attr::MOI.AnyAttribute) = MOI.supports(model.inner, attr)
MOI.supports_constraint(model::Optimizer, F::Type{<:MOI.VectorAffineFunction}, S::Type{<:MOI.AbstractVectorSet}) = MOI.supports_constraint(model.inner, F, S)


function _mult_subs(model, xterm::MOI.ScalarAffineTerm, ineq::MOI.ScalarAffineFunction)
    result = (xterm.coefficient * ineq.constant) * xterm.variable
    for yterm in ineq.terms
        push!(result.terms, MOI.ScalarAffineTerm(
            xterm.coefficient * yterm.coefficient,
            model.z[(xterm.variable, yterm.variable)],
        ))
    end
    return result
end
_mult_subs(model, term::MOI.ScalarAffineTerm, ineq::MOI.VariableIndex) = term.coefficient * model.z[(term.variable, ineq)]
_mult_subs(model::Optimizer, variableIndex::MOI.VariableIndex, ineq) = model.z[(variableIndex, ineq)]
function _mult_subs(model::Optimizer, affine::MOI.VectorAffineFunction{T}, ineq) where {T}
    result = MOI.Utilities.zero_with_output_dimension(MOI.VectorAffineFunction{T}, MOI.output_dimension(affine))
    for term in affine.terms
        scalar = _mult_subs(model, term.scalar_term, ineq)
        MOI.Utilities.operate_output_index!(+, T, term.output_index, result, scalar)
    end
    for i in eachindex(affine.constants)
        α = affine.constants[i]
        MOI.Utilities.operate_output_index!(+, T, i, result, α * ineq)
    end
    return result
end

# Check expressions contains only variables of the vertices or edges they are associated to
_check(_, e, _, ::MOI.VariableIndex, ::Nothing) = error("Expression `$e` is not assigned to any vertex or edge. You should set its `ConstraintVertexOrEdge` attribute if it is a constraint, or its `ObjectiveVertexOrEdge(*)` attribute if it is an objective function, where * is the vertex or edge identification.")
function _check(model, e, func, vi::MOI.VariableIndex, vertex::Int)
    variable_vertex = model.variable_vertex_or_edge[vi]
    if variable_vertex != vertex
        error("In expression `$e` of the vertex `$vertex`, the variable `$vi` of the function `$func` belongs to a different vertex `$variable_vertex`.")
    end
end
function _check(model, e, func, vi::MOI.VariableIndex, edge::Tuple{Int,Int})
    variable_vertex_or_edge = model.variable_vertex_or_edge[vi]
    if variable_vertex_or_edge != edge[1] && variable_vertex_or_edge != edge[2] && variable_vertex_or_edge != Graphs.Edge(edge...)
        error("In expression `$e` of the edge `$Graphs.Edge(edge...)`, the variable `$vi` of the function `$func` belongs to vertex or edge `$variable_vertex_or_edge` which is neither the source nor destination of the edge, nor the edge.")
    end
end
function _check(model, e, func::MOI.VectorAffineFunction, vertex_or_edge)
    for term in func.terms
        _check(model, e, func, term.scalar_term.variable, vertex_or_edge)
    end
end
function _check(model, e, vi::MOI.VariableIndex, vertex::Int)
    variable_vertex_or_edge = model.variable_vertex_or_edge[vi]
    if variable_vertex_or_edge != vertex
        error("In expression `$e` of the vertex `$vertex`, the variable `$vi` belongs to a different vertex or edge `$variable_vertex`.")
    end
end
function _check(model, e, vi::MOI.VariableIndex, edge::Tuple{Int,Int})
    variable_vertex_or_edge = model.variable_vertex_or_edge[vi]
    if variable_vertex_or_edge != edge[1] && variable_vertex_or_edge != edge[2] && variable_vertex_or_edge != Graphs.Edge(edge...)
        error("In expression `$e` of the edge `$Graphs.Edge(edge...)`, the variable `$vi` belongs to vertex or edge `$variable_vertex_or_edge` which is neither the source nor destination of the edge, nor the edge.")
    end
end

filter_attributes(attr) = true
filter_attributes(attr::VariableVertexOrEdge) = false
filter_attributes(attr::ConstraintVertexOrEdge) = false
filter_attributes(attr::Problem) = false
filter_attributes(attr::ObjectiveVertexOrEdge) = false
filter_attributes(attr::MOI.ObjectiveFunction) = false

# Function barrier to work around the type instability when getting `F` and `S`
function _add_constraints(dest::Optimizer, src::MOI.ModelLike, index_map, ::Type{F}, ::Type{S}) where {F,S}
    cis = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
    for ci in cis
        func = MOI.Utilities.map_indices(index_map, MOI.get(src, MOI.ConstraintFunction(), ci))
        set = MOI.get(src, MOI.ConstraintSet(), ci)
        vertex_or_edge = MOI.get(src, ConstraintVertexOrEdge(), ci)
        _check(dest, ci, func, vertex_or_edge)
        if vertex_or_edge isa Tuple{Int,Int}
            vertex_or_edge = Graphs.Edge(vertex_or_edge...)
        end
        # left of (5.8c) or (5.8e)
        linear = _mult_subs(dest, func, dest.y[vertex_or_edge])
        MOI.add_constraint(dest.inner, linear, set)
        if vertex_or_edge isa Int
            # It is a vertex so we can use Lemma 5.1 of the thesis
            # right of (5.8c)
            v = vertex_or_edge
            linear = _mult_subs(dest, func, 1.0 - dest.y[v])
            MOI.add_constraint(dest.inner, linear, set)
            # (5.8d)
            for u in Graphs.inneighbors(dest.graph, v)
                linear = _mult_subs(dest, func, dest.y[Graphs.Edge(u, v)])
                MOI.add_constraint(dest.inner, linear, set)
            end
            for u in Graphs.outneighbors(dest.graph, v)
                linear = _mult_subs(dest, func, 1.0 - dest.y[Graphs.Edge(v, u)])
                MOI.add_constraint(dest.inner, linear, set)
            end
        end
    end
    filtered = MOI.Utilities.ModelFilter(filter_attributes, src)
    MOI.Utilities.pass_attributes(dest.inner, filtered, index_map, cis)
end


function _set_objective(dest::Optimizer, src::MOI.ModelLike, index_map)
    result = MOI.Utilities.zero(MOI.ScalarAffineFunction{Float64})
    for vertex_or_edge in [collect(Graphs.vertices(dest.graph)); collect(Graphs.edges(dest.graph))]
        if (vertex_or_edge isa Int) v_e = vertex_or_edge else v_e = (Graphs.src(vertex_or_edge), Graphs.dst(vertex_or_edge)) end
        func = MOI.Utilities.map_indices(index_map, MOI.get(src, ObjectiveVertexOrEdge(v_e)))
        if !isnothing(func)
            _check(dest, ObjectiveVertexOrEdge(v_e), func, v_e)
            MOI.Utilities.operate!(+, Float64, result, _mult_subs(dest, func, dest.y[vertex_or_edge]))
        end
    end
    MOI.set(dest.inner, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), result)
end



struct ShortestPathProblem
    source::Int
    target::Int
end

MOI.Utilities.map_indices(::Function, p::ShortestPathProblem) = p

function _build_graph(graph::Graphs.DiGraph, src::MOI.ModelLike, ::Type{F}, ::Type{S}) where {F,S}
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        vertex_or_edge = MOI.get(src, ConstraintVertexOrEdge(), ci)
        if vertex_or_edge isa Tuple{Int,Int} # It's an edge
            u, v = vertex_or_edge
            Graphs.add_edge!(graph, u, v)
        end
    end
end

function _constrain_admissible_subgraphs(model::Optimizer, spp::ShortestPathProblem)
    for v in Graphs.vertices(model.graph)
        y = model.y[v]
        inv = MOI.VariableIndex[model.y[Graphs.Edge(u, v)] for u in Graphs.inneighbors(model.graph, v)]
        MOI.add_constraint(
            model.inner,
            y - dot(ones(length(inv)), inv),
            MOI.EqualTo(float(v == spp.source)),
        )
        outv = MOI.VariableIndex[model.y[Graphs.Edge(v, u)] for u in Graphs.outneighbors(model.graph, v)]
        MOI.add_constraint(
            model.inner,
            y - dot(ones(length(outv)), outv),
            MOI.EqualTo(float(v == spp.target)),
        )
    end
    for x in model.x
        v_e = model.variable_vertex_or_edge[x]
        if v_e isa Int # It is a vertex
            z = model.z[(x, model.y[v_e])]
            inv = MOI.VariableIndex[model.z[(x, model.y[Graphs.Edge(u, v_e)])] for u in Graphs.inneighbors(model.graph, v_e)]
            MOI.add_constraint(
                model.inner,
                z - dot(ones(length(inv)), inv) - float(v_e == spp.source) * x,
                MOI.EqualTo(0.0),
            )
            outv = MOI.VariableIndex[model.z[(x, model.y[Graphs.Edge(v_e, u)])] for u in Graphs.outneighbors(model.graph, v_e)]
            MOI.add_constraint(
                model.inner,
                z - dot(ones(length(outv)), outv) - float(v_e == spp.target) * x,
                MOI.EqualTo(0.0),
            )
        end
    end
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    MOI.empty!(dest)
    vis = MOI.get(src, MOI.ListOfVariableIndices())
    dest.x = MOI.add_variables(dest.inner, length(vis))
    index_map = MOI.Utilities.IndexMap()
    for (vi, x) in zip(vis, dest.x)
        index_map[vi] = x
        v_e = MOI.get.(src, VariableVertexOrEdge(), vi)
        if v_e isa Tuple{Int, Int} # It is an edge
            dest.variable_vertex_or_edge[x] = Graphs.Edge(v_e...)
        else # It is a vertex
            dest.variable_vertex_or_edge[x] = v_e
        end
    end
    vertices = unique!(sort(filter(x -> x isa Int, collect(values(dest.variable_vertex_or_edge)))))
    @assert vertices[1] == 1
    @assert vertices[end] == length(vertices)
    dest.graph = Graphs.DiGraph(length(vertices))
    # We build the graph first so that we have all the constraints to be able to do (5.8d)
    # That also allows to create all `y` and `z` variables
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        _build_graph(dest.graph, src, F, S)
    end
    # Create `y` variables
    for vertex in Graphs.vertices(dest.graph)
        # The constrant index is also returned but we don't need it.
        # We can construct it like this if needed:
        # `MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(model.z[key].value)`
        dest.y[vertex], _ = MOI.add_constrained_variable(dest.inner, MOI.ZeroOne())
    end
    for edge in Graphs.edges(dest.graph)
        dest.y[edge], _ = MOI.add_constrained_variable(dest.inner, MOI.ZeroOne())
    end
    # Create `z` variables
    for x in dest.x
        v_e = dest.variable_vertex_or_edge[x]
        dest.z[(x, dest.y[v_e])] = MOI.add_variable(dest.inner)
        if v_e isa Int # It is a vertex
            for u in Graphs.inneighbors(dest.graph, v_e)
                dest.z[(x, dest.y[Graphs.Edge(u, v_e)])] = MOI.add_variable(dest.inner)
            end
            for u in Graphs.outneighbors(dest.graph, v_e)
                dest.z[(x, dest.y[Graphs.Edge(v_e, u)])] = MOI.add_variable(dest.inner)
            end
        end
    end
    filtered = MOI.Utilities.ModelFilter(filter_attributes, src)
    MOI.Utilities.pass_attributes(dest.inner, filtered, index_map, vis)
    MOI.Utilities.pass_attributes(dest.inner, filtered, index_map)
    _set_objective(dest, src, index_map)
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        _add_constraints(dest, src, index_map, F, S)
    end
    _constrain_admissible_subgraphs(dest, MOI.get(src, Problem()))
    return index_map
end

MOI.get(model::Optimizer, attr::MOI.AbstractModelAttribute) = MOI.get(model.inner, attr)
MOI.get(model::Optimizer, attr::MOI.AbstractVariableAttribute, v) = MOI.get(model.inner, attr, v)
MOI.get(model::Optimizer, attr::MOI.AbstractConstraintAttribute, c) = MOI.get(model.inner, attr, c)

struct SubGraph <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::SubGraph) = true

MOI.Utilities.map_indices(::MOI.Utilities.IndexMap, ::SubGraph, graph::Graphs.DiGraph) = graph

function MOI.get(model::Optimizer, ::SubGraph)
    graph = Graphs.DiGraph(Graphs.nv(model.graph))
    for edge in Graphs.edges(model.graph)
        if MOI.get(model.inner, MOI.VariablePrimal(), model.y[edge]) > 0.5
            Graphs.add_edge!(graph, edge)
        end
    end
    return graph
end