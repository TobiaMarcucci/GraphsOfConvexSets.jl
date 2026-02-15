include("../MOI_structures.jl")

import Graphs
import MathOptInterface as MOI

struct ShortestPathProblem
    source::Int
    target::Int
end

MOI.Utilities.map_indices(::Function, p::ShortestPathProblem) = p

function _constrain_admissible_subgraphs(model::Optimizer, spp::ShortestPathProblem)
    for v in Graphs.vertices(model.graph)
        y = model.y[v]
        inv = MOI.VariableIndex[model.y[(u, v)] for u in Graphs.inneighbors(model.graph, v)] # Get the y variables of the incoming edges to v
        MOI.add_constraint( # sum of incoming edges to v + (v == src) = y_v
            model.inner,
            y - dot(ones(length(inv)), inv),
            MOI.EqualTo(float(v == spp.source)),
        )
        outv = MOI.VariableIndex[model.y[(v, u)] for u in Graphs.outneighbors(model.graph, v)] # Get the y variables of the outgoing edges from v
        MOI.add_constraint( # sum of outgoing edges to v + (v == target) = y_v
            model.inner,
            y - dot(ones(length(outv)), outv),
            MOI.EqualTo(float(v == spp.target)),
        )
    end
    MOI.add_constraint(model.inner, model.y[spp.source], MOI.EqualTo(1.)) # y_s = 1
    MOI.add_constraint(model.inner, model.y[spp.target], MOI.EqualTo(1.)) # y_t = 1
    for z in model.z
        v_e = model.variable_vertex_or_edge[z]
        if v_e isa Int # It is a vertex
            inv = MOI.VariableIndex[model.z_v_e[(z, model.y[(u, v_e)])] for u in Graphs.inneighbors(model.graph, v_e)] # Get the z_e variables of the incoming edges to v_e
            MOI.add_constraint( # z = sum of z_e of incoming edges to v + z * (v_e == src)
                model.inner,
                z - dot(ones(length(inv)), inv) - float(v_e == spp.source) * z,
                MOI.EqualTo(0.0),
            )
            outv = MOI.VariableIndex[model.z_v_e[(z, model.y[(v_e, u)])] for u in Graphs.outneighbors(model.graph, v_e)] # Get the z_e variables of the outgoing edges from v_e
            MOI.add_constraint( # z = sum of z_e of outgoing edges from v + z * (v_e == target)
                model.inner,
                z - dot(ones(length(outv)), outv) - float(v_e == spp.target) * z,
                MOI.EqualTo(0.0),
            )
        end
    end
end