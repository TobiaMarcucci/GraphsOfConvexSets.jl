function MOI.set(
    model::MOI.ModelLike,
    attr::VertexOrEdge,
    bridge::MOI.Bridges.Constraint.VectorizeBridge,
    value,
)
    return MOI.set(
        model,
        attr,
        bridge.vector_constraint,
        value,
    )
end
