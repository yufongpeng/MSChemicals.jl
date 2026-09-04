using MassSpecChemicals
using Statistics, StatsBase, TypedTables, DataPipes, Dictionaries
const MSC = MassSpecChemicals
import MassSpecChemicals: 
    chemicalname, 
    chemicalelements, 
    chemicalformula, 
    chemicalabbr, 
    chemicalsmiles, 
    completescheme, 
    elementalscheme, 
    adductionscheme, 
    detectedchemical, 
    charge
import Base: copy, hash, == 
using Test

test_show(x) = show(IOBuffer(), x)
macro test_noerror(x)
    return quote
        try 
            $x
            true
        catch e
            false
        end
    end
end

include(joinpath("objects", "generic.jl"))
include(joinpath("objects", "customized.jl"))
include(joinpath("objects", "isotopologues.jl"))
include(joinpath("objects", "lcms.jl"))

@info "Running tests"

@testset "MassSpecChemicals.jl" begin
    include("basic_function.jl")
    include("attr.jl")
    include("isotopologues.jl")
    include("lcms.jl")
    include("utils.jl")
end
