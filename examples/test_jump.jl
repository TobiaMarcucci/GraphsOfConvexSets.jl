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

println(typeof(model) <: JuMP.AbstractModel, "\n")

v1 = Vertex(model, "pirlou")
v2 = Vertex(model, "youpidou")
e1 = Edge(model, (v1,v2))

@variable(v1, x)