using JuMP

const NameType = Union{Int,String}

struct Vertex{M <: JuMP.AbstractModel} <: JuMP.AbstractModel
    model::M
    vertex::NameType
end

struct Edge{M <: JuMP.AbstractModel} <: JuMP.AbstractModel
    model::M
    edge::Tuple{Vertex, Vertex}
end

JuMP.object_dictionary(v::Vertex) = JuMP.object_dictionary(v.model)
JuMP.object_dictionary(e::Edge) = JuMP.object_dictionary(e.model)