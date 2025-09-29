module GraphsOfConvexSets

include("MOI_wrapper.jl")
include("bridges.jl")
include("graph_problems/problem_structures.jl")

# Taken from JuMP.jl, exports all symbols not starting with `_`
const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__; all = true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS ||
       startswith(sym_string, "_") ||
       startswith(sym_string, "@_")
        continue
    end
    if !(
        Base.isidentifier(sym) ||
        (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end]))
    )
        continue
    end
    @eval export $sym
end

end
