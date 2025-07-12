import Graphs
import MathOptInterface as MOI
using LinearAlgebra

const VertexOrEdgeType = Union{Int,Graphs.Edge{Int}}

struct Problem <: MOI.AbstractModelAttribute end
struct VariableVertex <: MOI.AbstractVariableAttribute end
struct VertexOrEdge <: MOI.AbstractConstraintAttribute end

mutable struct Optimizer{O} <: MOI.AbstractOptimizer
    inner::O
    graph::Union{Nothing,Graphs.SimpleDiGraph{Int}}
    variable_vertex::Union{Nothing,Vector{Int}}
    x::Union{Nothing,Vector{MOI.VariableIndex}}
    y::Dict{VertexOrEdgeType,MOI.VariableIndex}
    z::Dict{Tuple{MOI.VariableIndex,MOI.VariableIndex},MOI.VariableIndex} # (x, y) -> z
    #vertex_or_edge::Dict{MOI.ConstraintIndex,VertexOrEdgeType}
    function Optimizer(optimizer_constructor)
        inner = MOI.instantiate(optimizer_constructor, with_bridge_type=Float64)
        return new{typeof(inner)}(
            inner,
            nothing,
            nothing,
            nothing,
            Dict{VertexOrEdgeType,MOI.VariableIndex}(),
            Dict{Tuple{MOI.VariableIndex,MOI.VariableIndex},MOI.VariableIndex}(),
            #Dict{MOI.ConstraintIndex,VertexOrEdgeType}(),
        )
    end
end

MOI.is_empty(model::Optimizer) = MOI.is_empty(model.inner)
function MOI.empty!(model::Optimizer)
    MOI.empty!(model.inner)
    model.graph = nothing
    model.variable_vertex = nothing
    model.x = nothing
    empty!(model.y)
    empty!(model.z)
    #empty!(model.vertex_or_edge)
    return
end
MOI.optimize!(model::Optimizer) = MOI.optimize!(model.inner)
function MOI.supports(model::Optimizer, attr::MOI.AnyAttribute)
    return MOI.supports(model.inner, attr)
end
function MOI.supports_constraint(model::Optimizer, F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})
    return MOI.supports_constraint(model.inner, F, S)
end

#function _homogenize(model::Optimizer, constant::Number, vertex_or_edge)
#    key = vertex_or_edge
#    if haskey(model.y, key)
#        # The constrant index is also returned but we don't need it.
#        # We can construct it like this if needed:
#        # `MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(model.z[key].value)`
#        model.y[key], _ = MOI.add_constrained_variable(model, MOI.ZeroOne())
#    end
#    return MOI.ScalarAffineTerm(constant, model.y[key])
#end
#
#function _homogenize(model::Optimizer, constants::Vector, vertex_or_edge)
#    return MOI.VectorAffineTerm.(eachindex(constants), _homogenize.(model, constants, vertex_or_edge))
#end
#
#function _homogenize(model::Optimizer, vi::MOI.VariableIndex, vertex_or_edge)
#    key = (vertex_or_edge, vi)
#    if haskey(model.z, key)
#        # The constrant index is also returned but we don't need it.
#        # We can construct it like this if needed:
#        # `MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(model.z[key].value)`
#        model.z[key] = MOI.add_variable(model)
#    end
#    return model.z[key]
#end
#
#function _homogenize(model::Optimizer, term::MOI.ScalarAffineTerm, vertex_or_edge)
#    return MOI.ScalarAffineTerm(
#        term.coefficient,
#        _homogenize(model, term.variable, vertex_or_edge),
#    )
#end

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

function _mult_subs(model, term::MOI.ScalarAffineTerm, ineq::MOI.VariableIndex)
    return term.coefficient * model.z[(term.variable, ineq)]
end

# Return the result of `affine * ineq` after substituting the bilinear terms `x * y` into `z`
function _mult_subs(model::Optimizer, affine::MOI.ScalarAffineFunction{T}, ineq) where {T}
    result = zero(MOI.ScalarAffineFunction{T})
    for term in affine.terms
        scalar = _mult_subs(model, term, ineq)
        MOI.Utilities.operate!(+, T, result, scalar)
    end
    MOI.Utilities.operate!(+, T, result, affine.constant * ineq)
    return result
end

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

function _build_graph(graph::Graphs.DiGraph, src::MOI.ModelLike, ::Type{F}, ::Type{S}) where {F,S}
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        vertex_or_edge = MOI.get(src, VertexOrEdge(), ci)
        if vertex_or_edge isa Tuple{Int,Int} # It's an edge
            u, v = vertex_or_edge
            Graphs.add_edge!(graph, u, v)
        end
    end
end

# Function barrier to work around the type instability when getting `F` and `S`
function _add_constraints(dest::Optimizer, src::MOI.ModelLike, ::Type{F}, ::Type{S}) where {F,S}
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        func = MOI.get(src, MOI.ConstraintFunction(), ci)
        set = MOI.get(src, MOI.ConstraintSet(), ci)
        @assert !(set isa MOI.Interval) # TODO exclude with `supports_constraint`
        if MOI.Utilities.supports_shift_constant(typeof(set))
            # We need to move the constant to the set to homogenize it
            constant = MOI.constant(set)
            func -= constant
            set = MOI.Utilities.shift_constant(set, -constant)
        end
        vertex_or_edge = MOI.get(src, VertexOrEdge(), ci)
        # left of (5.8c) or (5.8e)
        linear = _mult_subs(dest, func, dest.y[vertex_or_edge])
        _add_constraint(dest.inner, linear, set)
        if vertex_or_edge isa Int
            # It is a vertex so we can use Lemma 5.1 of the thesis
            # right of (5.8c)
            v = vertex_or_edge
            linear = _mult_subs(dest, func, 1.0 - dest.y[v])
            _add_constraint(dest.inner, linear, set)
            # (5.8d)
            for u in Graphs.inneighbors(dest.graph, v)
                linear = _mult_subs(dest, func, dest.y[(u, v)])
                _add_constraint(dest.inner, linear, set)
            end
            for u in Graphs.outneighbors(dest.graph, v)
                linear = _mult_subs(dest, func, 1.0 - dest.y[(v, u)])
                _add_constraint(dest.inner, linear, set)
            end
        end
    end
end

# We need to put the constant back in the set when appropriate
_add_constraint(model, f::MOI.AbstractScalarFunction, s) = MOI.Utilities.normalize_and_add_constraint(model, f, s)
_add_constraint(model, f::MOI.AbstractVectorFunction, s) = MOI.add_constraint(model, f, s)

struct ShortestPathProblem
    source::Int
    target::Int
end

MOI.Utilities.map_indices(::Function, p::ShortestPathProblem) = p

function _constrain_admissible_subgraphs(model::Optimizer, spp::ShortestPathProblem)
    for v in Graphs.vertices(model.graph)
        y = model.y[v]
        inv = MOI.VariableIndex[model.y[(u, v)] for u in Graphs.inneighbors(model.graph, v)]
        MOI.add_constraint(
            model.inner,
            y - dot(ones(length(inv)), inv),
            MOI.EqualTo(float(v == spp.source)),
        )
        outv = MOI.VariableIndex[model.y[(v, u)] for u in Graphs.outneighbors(model.graph, v)]
        MOI.add_constraint(
            model.inner,
            y - dot(ones(length(outv)), outv),
            MOI.EqualTo(float(v == spp.target)),
        )
    end
    for i in eachindex(model.x)
        x = model.x[i]
        v = model.variable_vertex[i]
        z = model.z[(x, model.y[v])]
        inv = MOI.VariableIndex[model.z[(x, model.y[Graphs.Edge(u, v)])] for u in Graphs.inneighbors(model.graph, v)]
        MOI.add_constraint(
            model.inner,
            z - dot(ones(length(inv)), inv) - float(v == spp.source) * x,
            MOI.EqualTo(0.0),
        )
        outv = MOI.VariableIndex[model.z[(x, model.y[Graphs.Edge(v, u)])] for u in Graphs.outneighbors(model.graph, v)]
        MOI.add_constraint(
            model.inner,
            z - dot(ones(length(outv)), outv) - float(v == spp.target) * x,
            MOI.EqualTo(0.0),
        )
    end
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    MOI.empty!(dest)
    vis = MOI.get(src, MOI.ListOfVariableIndices())
    dest.x = MOI.add_variables(dest.inner, length(vis))
    index_map = MOI.Utilities.IndexMap()
    for i in eachindex(vis)
        index_map[vis[i]] = dest.x[i]
    end
    dest.variable_vertex = MOI.get.(src, VariableVertex(), vis)
    vertices = unique!(sort(dest.variable_vertex))
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
    for i in eachindex(dest.x)
        x = dest.x[i]
        v = dest.variable_vertex[i]
        dest.z[(x, dest.y[v])] = MOI.add_variable(dest.inner)
        for u in Graphs.inneighbors(dest.graph, v)
            dest.z[(x, dest.y[Graphs.edge(u, v)])] = MOI.add_variable(dest.inner)
        end
        for u in Graphs.outneighbors(dest.graph, v)
            dest.z[(x, dest.y[Graphs.edge(v, u)])] = MOI.add_variable(dest.inner)
        end
    end
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        _add_constraints(dest, src, F, S)
    end
    _constrain_admissible_subgraphs(dest, MOI.get(src, Problem()))
    return index_map
end
