import Graphs
import MathOptInterface as MOI

const VertexOrEdgeType = Union{Int,Tuple{Int,Int}}

# Variables, Constraints and Model attributes
struct VariableVertexOrEdge <: MOI.AbstractVariableAttribute end     # Associate to a variable its vertex or edge
struct ConstraintVertexOrEdge <: MOI.AbstractConstraintAttribute end # Associate to a constraint its vertex or edge
struct ObjectiveVertexOrEdge <: MOI.AbstractModelAttribute           # Associate to a vertex or edge an objective
    vertex_or_edge::VertexOrEdgeType
end
struct Problem <: MOI.AbstractModelAttribute end                     # Associate a problem to the model, e.g. SPP

# Optimizer
mutable struct Optimizer{O} <: MOI.AbstractOptimizer
    inner::O # Optimizer to use to solve the optimization problem
    graph::Union{Nothing,Graphs.SimpleDiGraph{Int}} # Represents the graph
    variable_vertex_or_edge::Dict{MOI.VariableIndex,VertexOrEdgeType} # associate to each variable z a vertex or edge
    z::Union{Nothing, Vector{MOI.VariableIndex}} # Variables z_v and z_e of the model
    y::Dict{VertexOrEdgeType,MOI.VariableIndex}  # Associate to each vertex or edge a binary variable y
    z_v_e::Dict{Tuple{MOI.VariableIndex, MOI.VariableIndex}, MOI.VariableIndex} # associate to each variable z_v and variable y_e, with e adjacent to v, a variable z_v_e : (z, y) => z_v_e
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

MOI.is_empty(model::Optimizer) = MOI.is_empty(model.inner) # Returns if the model is empty
function MOI.empty!(model::Optimizer) # Empty the model
    MOI.empty!(model.inner)
    model.graph = nothing
    empty!(model.variable_vertex_or_edge)
    model.z = nothing
    empty!(model.y)
    empty!(model.z_v_e)
    return
end

MOI.supports(model::Optimizer, attr::MOI.AnyAttribute) = MOI.supports(model.inner, attr)
MOI.supports_constraint(model::Optimizer, F::Type{<:MOI.VectorAffineFunction}, S::Type{<:MOI.AbstractVectorSet}) = MOI.supports_constraint(model.inner, F, S)

MOI.get(model::Optimizer, attr::MOI.AbstractModelAttribute) = MOI.get(model.inner, attr)
MOI.get(model::Optimizer, attr::MOI.AbstractVariableAttribute, v) = MOI.get(model.inner, attr, v)
MOI.get(model::Optimizer, attr::MOI.AbstractConstraintAttribute, c) = MOI.get(model.inner, attr, c)

MOI.optimize!(model::Optimizer) = MOI.optimize!(model.inner)