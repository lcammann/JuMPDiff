module JuMPDiff

using Reexport;
    @reexport using JuMP

export model_jac, model_jac!, model_hess;

include("Functions.jl")

end # module
