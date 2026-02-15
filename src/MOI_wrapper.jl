using LinearAlgebra
using JuMP
import Graphs
import MathOptInterface as MOI

include("graph_problems/shortest_path.jl")

###################################################################

const VertexOrEdgeType = Union{Int,Tuple{Int,Int}}

# Variables, Constraints and Model attributes
struct VariableVertexOrEdge <: MOI.AbstractVariableAttribute end     # Associate to a variable its vertex or edge
struct ConstraintVertexOrEdge <: MOI.AbstractConstraintAttribute end # Associate to a constraint its vertex or edge
struct VertexOrEdgeObjective <: MOI.AbstractModelAttribute           # Associate to a vertex or edge an objective
    vertex_or_edge::VertexOrEdgeType
end
struct Problem <: MOI.AbstractModelAttribute end                     # Associate a problem to the model, e.g. SPP

# Optimizer
mutable struct GCSOptimizer{O} <: MOI.AbstractOptimizer
    inner::O # Optimizer to use to solve the optimization problem
    graph::Union{Nothing,Graphs.SimpleDiGraph{Int}} # Represents the graph
    variable_vertex_or_edge::Dict{MOI.VariableIndex,VertexOrEdgeType} # associate to each variable z a vertex or edge
    z::Union{Nothing, Vector{MOI.VariableIndex}} # Variables z_v and z_e of the model
    y::Dict{VertexOrEdgeType,MOI.VariableIndex}  # Associate to each vertex or edge a binary variable y
    z_v_e::Dict{Tuple{MOI.VariableIndex, MOI.VariableIndex}, MOI.VariableIndex} # associate to each variable z_v and variable y_e, with e adjacent to v, a variable z_v_e : (z, y) => z_v_e
    function GCSOptimizer(optimizer_constructor)
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

###################################################################

"""
_mult_subs is a function that takes as input an affine and a linear functions:

    - one      of the kind   Az + b;

    - a second of the kind   a y_v + ∑_(e ∈ I_v) b_e y_e for some vertex v

and that returns their product, substituting the following bilinear terms:

    - z_v y_v = z_v

    - z_v y_e = z_v^e

    - z_e y_v = z_e y_e = z_e
"""
function _mult_subs(model::Optimizer, zterm::MOI.ScalarAffineTerm{T}, ineq::MOI.ScalarAffineFunction) where {T}
    # This method consider the case where the first affine function is scalar with one linear term of the kind a*z
    result = MOI.Utilities.zero(MOI.ScalarAffineFunction{T})
    for yterm in ineq.terms
        push!(result.terms, MOI.ScalarAffineTerm(
            zterm.coefficient * yterm.coefficient,
            haskey(model.z_v_e, (zterm.variable, yterm.variable)) ? model.z_v_e[(zterm.variable, yterm.variable)] : zterm.variable,
        ))
    end
    return result
end

# This method consider the case where the first affine function is scalar with one linear term of the kind a*z, and the second linear function is just y_v or y_e.
_mult_subs(model::Optimizer, term::MOI.ScalarAffineTerm, ineq::MOI.VariableIndex) = term.coefficient * (haskey(model.z_v_e, (term.variable, ineq)) ? model.z_v_e[(term.variable, ineq)] : term.variable)

function _mult_subs(model::Optimizer, variableIndex::MOI.VariableIndex, ineq::MOI.ScalarAffineFunction)
    # This method consider the case where the first affine function is just z_v or z_e
    result = MOI.Utilities.zero(MOI.ScalarAffineFunction{Float64})
    for yterm in ineq.terms
        push!(result.terms, MOI.ScalarAffineTerm(
            yterm.coefficient,
            haskey(model.z_v_e, (variableIndex, yterm.variable)) ? model.z_v_e[(variableIndex, yterm.variable)] : variableIndex))
    end
    return result
end

# This method consider the case where the first affine function is just z_v or z_e, and the second linear function is just y_v or y_e.
_mult_subs(model::Optimizer, variableIndex::MOI.VariableIndex, ineq::MOI.VariableIndex) = haskey(model.z_v_e, (variableIndex, ineq)) ? model.z_v_e[(variableIndex, ineq)] : variableIndex

function _mult_subs(model::Optimizer, affine::MOI.ScalarAffineFunction{T}, ineq) where {T}
    # This method handles a whole scalar affine function a*z + b, taking the terms one by one.
    result = MOI.Utilities.zero(MOI.ScalarAffineFunction{T})
    for term in affine.terms
        scalar = _mult_subs(model, term, ineq)
        MOI.Utilities.operate!(+, T, result, scalar)
    end
    MOI.Utilities.operate!(+, T, result, affine.constant * ineq) # product of the cst term b with the linear function ineq
    return result
end

function _mult_subs(model::Optimizer, affine::MOI.VectorAffineFunction{T}, ineq) where {T}
    # This method handles a vector affine function Az + b, taking the terms of each scalar affine expression one by one.
    result = MOI.Utilities.zero_with_output_dimension(MOI.VectorAffineFunction{T}, MOI.output_dimension(affine))
    for term in affine.terms
        scalar = _mult_subs(model, term.scalar_term, ineq)
        MOI.Utilities.operate_output_index!(+, T, term.output_index, result, scalar)
    end
    for i in eachindex(affine.constants) # product of the cst term b with the linear function ineq
        α = affine.constants[i]
        MOI.Utilities.operate_output_index!(+, T, i, result, α * ineq)
    end
    return result
end

########################################################################################################

"""
Check expressions contain only variables of the vertices or edges they are associated to
"""
_check(_, e, ::MOI.VariableIndex, ::Nothing) = error("Expression `$e` is not assigned to any vertex or edge. You should set its `ConstraintVertexOrEdge` attribute if it is a constraint, or its `VertexOrEdgeObjective(*)` attribute if it is an objective function, where * is the vertex or edge identification.")

function _check(model::Optimizer, e, vi::MOI.VariableIndex, vertex::Int)
    # Check if a variable is associated to a vertex
    variable_vertex_or_edge = model.variable_vertex_or_edge[vi]
    if variable_vertex_or_edge != vertex
        error("In expression `$e` of the vertex `$vertex`, the variable `$vi` belongs to a different vertex or edge `$variable_vertex`.")
    end
end

function _check(model::Optimizer, e, vi::MOI.VariableIndex, edge::Tuple{Int,Int})
    # Check if a variable is associated to an edge or its incident vertices.
    variable_vertex_or_edge = model.variable_vertex_or_edge[vi]
    if variable_vertex_or_edge != edge[1] && variable_vertex_or_edge != edge[2] && variable_vertex_or_edge != edge
        error("In expression `$e` of the edge `$Graphs.Edge(edge...)`, the variable `$vi` belongs to vertex or edge `$variable_vertex_or_edge` which is neither the source nor destination of the edge, nor the edge.")
    end
end

function _check(model, e, func::MOI.ScalarAffineFunction, vertex_or_edge)
    # Check for a scalar affine function by looking at the variables of all the terms separately
    for term in func.terms
        _check(model, e, term.variable, vertex_or_edge)
    end
end

function _check(model, e, func::MOI.VectorAffineFunction, vertex_or_edge)
    # Check for a vector affine function by looking at the variables of all the terms separately
    for term in func.terms
        _check(model, e, term.scalar_term.variable, vertex_or_edge)
    end
end

#############################################################################################

"""
function to transfer the variables, constraint and model attributes of the src model to the dest model, without the attributes specific to a GCS
"""
filter_attributes(attr) = true # Default method
filter_attributes(attr::VariableVertexOrEdge) = false
filter_attributes(attr::ConstraintVertexOrEdge) = false
filter_attributes(attr::Problem) = false
filter_attributes(attr::VertexOrEdgeObjective) = false
filter_attributes(attr::MOI.ObjectiveFunction) = false

###########################################################################################

# Function barrier to work around the type instability when getting `F` and `S`
function _add_constraints(dest::Optimizer, src::MOI.ModelLike, index_map, ::Type{F}, ::Type{S}) where {F,S}
    cis = MOI.get(src, MOI.ListOfConstraintIndices{F,S}()) # Get all the constraints of the type f in s, with f ∈ F a function and s ∈ S a set
    for ci in cis
        func = MOI.Utilities.map_indices(index_map, MOI.get(src, MOI.ConstraintFunction(), ci)) # get the function f, and map the x variables of the src model to the z variables of the dest model
        set  = MOI.get(src, MOI.ConstraintSet(), ci) # get the constraint set s
        vertex_or_edge = MOI.get(src, ConstraintVertexOrEdge(), ci) # get the vertex or edge associated to the constraint
        _check(dest, ci, func, vertex_or_edge) # check if all the variables of func belong to the vertex or edge (or vertices incident to the edge)

        # Recall that for a conic constraint of the kind Ax + b ∈ K, the homogenized constraint is A(z*y) + by ∈ K, substituting z*y by the corresponding variable
        if vertex_or_edge isa Tuple{Int,Int} # It is an edge, so we add the constraint A(z*y_e) + by_e ∈ K - corresponds to (14.c)
            linear = _mult_subs(dest, func, dest.y[vertex_or_edge]) # Computes A(z*y_e) + by_e
            MOI.add_constraint(dest.inner, linear, set)
        else # It is a vertex, so we add the constraints A(z_v*y_e) + b*y_e ∈ K and Az_v(y_v - y_e) + b(y_v - y_e) ∈ K - corresponds to (14.d) and (14.e)
            v = vertex_or_edge
            for u in Graphs.inneighbors(dest.graph, v)
                e = (u,v)
                # (14.d)
                linear = _mult_subs(dest, func, dest.y[e]) # Computes A(z_v*y_e) + b*y_e
                MOI.add_constraint(dest.inner, linear, set)
                # (14.e)
                linear = _mult_subs(dest, func, 1.0dest.y[v] - 1.0dest.y[e]) # Computes Az_v(y_v - y_e) + b(y_v - y_e)
                MOI.add_constraint(dest.inner, linear, set)
            end
            for u in Graphs.outneighbors(dest.graph, v)
                e = (v,u)
                # (14.d)
                linear = _mult_subs(dest, func, dest.y[e]) # Computes A(z_v*y_e) + b*y_e
                MOI.add_constraint(dest.inner, linear, set)
                # (14.e)
                linear = _mult_subs(dest, func, 1.0dest.y[v] - 1.0dest.y[e]) # Computes Az_v(y_v - y_e) + b(y_v - y_e)
                MOI.add_constraint(dest.inner, linear, set)
            end
        end
    end

    filtered = MOI.Utilities.ModelFilter(filter_attributes, src)
    MOI.Utilities.pass_attributes(dest.inner, filtered, index_map, cis) # Pass the attributes of the constraints from the src model to the dest model, according to the filter_attributes function
end

#################################################################################################

"""
Function used by the user to set the objective of a vertex or edge
"""
function set_vertex_or_edge_objective(model::Model, v_e::VertexOrEdgeType, func::Union{Number, AbstractVariableRef, GenericAffExpr})
    MOI.set(model, VertexOrEdgeObjective(v_e), moi_function(func))
end

"""
Construct the objective of the dest model, as the sum of the homogenized costs of the vertices and edges.
If the cost associated to a vertex or edge is of the form ax + b, its homogenized version is az + by
"""
function _set_objective(dest::Optimizer, src::MOI.ModelLike, index_map)
    result = MOI.Utilities.zero(MOI.ScalarAffineFunction{Float64})
    for vertex_or_edge in [collect(Graphs.vertices(dest.graph)); collect(Graphs.edges(dest.graph))]
        v_e = vertex_or_edge isa Int ? vertex_or_edge : (Graphs.src(vertex_or_edge), Graphs.dst(vertex_or_edge))
        func = MOI.Utilities.map_indices(index_map, MOI.get(src, VertexOrEdgeObjective(v_e)))
        if !isnothing(func)
            _check(dest, VertexOrEdgeObjective(v_e), func, v_e)
            MOI.Utilities.operate!(+, Float64, result, _mult_subs(dest, func, dest.y[v_e]))
        end
    end
    MOI.set(dest.inner, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), result)
end

##############################################################################################################

"""
Build the edges of the graph based on the constraints binding vertices
"""
function _build_graph(graph::Graphs.DiGraph, src::MOI.ModelLike, ::Type{F}, ::Type{S}) where {F,S}
    for ci in MOI.get(src, MOI.ListOfConstraintIndices{F,S}()) # Get all the constraints of the type f in s, with f ∈ F a function and s ∈ S a set
        vertex_or_edge = MOI.get(src, ConstraintVertexOrEdge(), ci)
        if vertex_or_edge isa Tuple{Int,Int} # The constraint is on an edge --> we add the edge to the graph
            u, v = vertex_or_edge
            Graphs.add_edge!(graph, u, v)
        end
    end
end

#################################################################################################################

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike)
    MOI.empty!(dest)

    # Create the z_v and z_e variables
    vis = MOI.get(src, MOI.ListOfVariableIndices()) # Get a list of the variable indices of the src model
    dest.z = MOI.add_variables(dest.inner, length(vis)) # Create a z variable for the dest model for each variable in vis
    index_map = MOI.Utilities.IndexMap() # Map between the x variables of the src model and the z variables of the dest model
    for (vi, z) in zip(vis, dest.z)
        index_map[vi] = z # Associate each variable of the src model to its corresponding variable in the dest model
        dest.variable_vertex_or_edge[z] = MOI.get.(src, VariableVertexOrEdge(), vi) # Associate each variable to its vertex or edge
    end

    # Construct the graph
    vertices = unique!(sort(filter(x -> x isa Int, collect(values(dest.variable_vertex_or_edge))))) # Get the vertices names, that should be 1 to n with n the number of vertices
    @assert vertices[1] == 1
    @assert vertices[end] == length(vertices)
    dest.graph = Graphs.DiGraph(length(vertices)) # Create a directed graph with n vertices
    # We build the graph first so that we have all the constraints to be able to do (5.8d)
    # That also allows to create all `y` and `z_v_e` variables
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        _build_graph(dest.graph, src, F, S) # Add the edges to the graph
    end

    # Create `y` variables
    for vertex in Graphs.vertices(dest.graph) # y_v variables
        # The constrant index is also returned but we don't need it.
        # We can construct it like this if needed:
        # `MOI.ConstraintIndex{MOI.VariableIndex,MOI.ZeroOne}(model.z[key].value)`
        dest.y[vertex], _ = MOI.add_constrained_variable(dest.inner, MOI.ZeroOne())
    end
    for edge in Graphs.edges(dest.graph) # y_e variables
        dest.y[(Graphs.src(edge), Graphs.dst(edge))], _ = MOI.add_constrained_variable(dest.inner, MOI.ZeroOne())
    end
    # Create `z_v_e` variables
    for z in dest.z
        v_e = dest.variable_vertex_or_edge[z]
        if v_e isa Int # It is a vertex
            for u in Graphs.inneighbors(dest.graph, v_e) # Incoming edges
                dest.z_v_e[(z, dest.y[(u, v_e)])] = MOI.add_variable(dest.inner)
            end
            for w in Graphs.outneighbors(dest.graph, v_e) # Outgoing edges
                dest.z_v_e[(z, dest.y[(v_e, w)])] = MOI.add_variable(dest.inner)
            end
        end
    end

    # Copy attributes of the src model to the dest model
    filtered = MOI.Utilities.ModelFilter(filter_attributes, src)
    MOI.Utilities.pass_attributes(dest.inner, filtered, index_map, vis)
    MOI.Utilities.pass_attributes(dest.inner, filtered, index_map)

    _set_objective(dest, src, index_map) # Set the objective of the dest model to the sum of the objectives on all chosen edges and vertices
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent()) # Set the constraints of the dest model
        _add_constraints(dest, src, index_map, F, S)
    end

    # Set the constraints specific to the problem
    _constrain_admissible_subgraphs(dest, MOI.get(src, Problem()))
    
    return index_map
end

####################################################################################################

struct SubGraph <: MOI.AbstractModelAttribute end
MOI.is_set_by_optimize(::SubGraph) = true

MOI.Utilities.map_indices(::MOI.Utilities.IndexMap, ::SubGraph, graph::Graphs.DiGraph) = graph

function MOI.get(model::Optimizer, ::SubGraph)
    graph = Graphs.DiGraph(Graphs.nv(model.graph))
    for edge in Graphs.edges(model.graph)
        edge = (Graphs.src(edge), Graphs.dst(edge))
        if MOI.get(model.inner, MOI.VariablePrimal(), model.y[edge]) > 0.5
            Graphs.add_edge!(graph, edge)
        end
    end
    return graph
end

####################################################################################################

function JuMP.add_variable(v::Vertex, var...)
    var_ref = JuMP.add_variable(v.model, var...)
    MOI.set(v.model, VariableVertexOrEdge(), var_ref, v.vertex)
end

function JuMP.add_constraint(v::Vertex, con)
    con_ref = JuMP.add_constraint(v.model, con)
    MOI.set(v.model, ConstraintVertexOrEdge(), con_ref, v.vertex)
end

function JuMP.add_variable(e::Edge, var)
    var_ref = JuMP.add_variable(e.model, var)
    MOI.set(e.model, VariableVertexOrEdge(), var_ref, e.edge)
end

function JuMP.add_constraint(e::Edge, con)
    con_ref = JuMP.add_constraint(e.model, con)
    MOI.set(e.model, ConstraintVertexOrEdge(), con_ref, e.edge)
end