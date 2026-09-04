"""
    Chemical <: AbstractChemical

Unstructured chemical type with its name, elements (formula), and additional properties.

# Fields 
* `name::String`: a unique chemical name.
* `elements::Vector{Pair{String, Int}}`: chemical elements.
* `property::Vector{Pair{Symbol, Any}}`: additional properties; the pairs repressent names and values.

# Constructors
    Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}, property::Vector{Pair{Symbol, Any}})
    Chemical(name::AbstractString, elements::Dict{String, Int}, property::Vector{Pair{Symbol, Any}})
    Chemical(name::AbstractString, formula::AbstractString, property::Vector{Pair{Symbol, Any}})
    Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}; kwargs...)
    Chemical(name::AbstractString, elements::Dict{String, Int}; kwargs...)
    Chemical(name::AbstractString, formula::AbstractString; kwargs...)

`kwargs` are collected into field `property`.
"""
struct Chemical <: AbstractChemical
    name::String
    elements::Vector{Pair{String, Int}}
    property::Vector{Pair{Symbol, Any}}
end

Chemical(name::AbstractString, elements; kwargs...) = Chemical(name, elements, collect(kwargs))
Chemical(name::AbstractString, formula::AbstractString, property) = Chemical(name, chemicalelements(formula), property)
Chemical(name::AbstractChemical, elements::Dict, property) = Chemical(name, collect(elements), property)

"""
    FormulaChemical <: AbstractChemical

Unstructured chemical type with elements (formula), and additional properties. Chemical name will be formula.

# Fields 
* `elements::Vector{Pair{String, Int}}`: chemical elements.
* `property::Vector{Pair{Symbol, Any}}`: additional properties; the pairs repressent names and values.

# Constructors
    FormulaChemical(elements::Vector{Pair{String, Int}}, property::Vector{Pair{Symbol, Any}})
    FormulaChemical(elements::Dict{String, Int}, property::Vector{Pair{Symbol, Any}})
    FormulaChemical(formula::AbstractString, property::Vector{Pair{Symbol, Any}})
    FormulaChemical(elements::Vector{Pair{String, Int}}; kwargs...)
    FormulaChemical(elements::Dict{String, Int}; kwargs...)
    FormulaChemical(formula::AbstractString; kwargs...)

`kwargs` are collected into field `property`.
"""
struct FormulaChemical <: AbstractChemical
    elements::Vector{Pair{String, Int}}
    property::Vector{Pair{Symbol, Any}}
end

FormulaChemical(elements; kwargs...) = FormulaChemical(elements, collect(kwargs))
FormulaChemical(formula::AbstractString, property) = FormulaChemical(chemicalelements(formula), property)
FormulaChemical(elements::Dict, property) = FormulaChemical(collect(elements), property)

"""
    ChemicalTransition{T<:AbstractChemicalsSchema} <: AbstractChemical

Chemical transition in MSⁿ. Products can be any subtype of `AbstractChemicalsSchema` representing chemical entity or species.

# Fields 
* `transition::Vector{T}`.

# Constructors
    ChemicalTransition(transition::Vector)
    ChemicalTransition(precursor, products...)`
    
`products` are pushed into `precursor` to construct `transition`.
"""
struct ChemicalTransition{T<:AbstractChemicalsSchema} <: AbstractChemical
    transition::Vector{T}
end

function ChemicalTransition(ct...) 
    ChemicalTransition(mapreduce(_transition, vcat, ct))
end

_transition(x::ChemicalTransition) = chemicaltransition(x)
_transition(x) = x

"""
    Isobars{T<:AbstractChemical, N} <: AbstractChemical

Chemicals with similar m/z.

# Fields 
* `chemicals::Vector{T}`: a vector of chemicals.
* `abundnace::VecOrMat{N}`: the abundance of each chemical. If chemicals are trasitions, this should be a matrix, and each column is the abundance of each ms stage.

# Constructors 
    Isobars(chemicals::Vector, abundance::Vector)
    Isobars(chemicals::Vector{<:ChemicalTransition}, abundance::Vector)
    Isobars(chemicals::Vector{<:ChemicalTransition}, abundance::Matrix)
"""
struct Isobars{T<:AbstractChemical, N} <: AbstractChemical
    chemicals::Vector{T}
    abundance::VecOrMat{N}
    function Isobars(chemicals::Vector{T}, abundance::Vector{N}) where {T, N}
        id = sortperm(abundance; rev = true)
        new{T, N}(chemicals[id], abundance[id])
    end
    function Isobars(chemicals::Vector{T}, abundance::Vector{N}) where {T<:ChemicalTransition, N}
        allequal(msstage, chemicals) || throw(ArgumentError("All chemicals should have the same `msstage`."))
        length(chemicals) == length(abundance) || throw(ArgumentError("The length of `chemicals` and `abundance` are not equaled."))
        id = sortperm(abundance; rev = true)
        ab = hcat([abundance[id] for _ in 1:msstage(chemicals[begin])]...)
        new{T, N}(chemicals[id], ab)
    end
    function Isobars(chemicals::Vector{T}, abundance::Vector{Vector{N}}) where {T<:ChemicalTransition, N}
        allequal(msstage, chemicals) || throw(ArgumentError("All chemicals should have the same `msstage`."))
        length(abundance) == msstage(chemicals[begin]) || throw(ArgumentError("The length of `abundance` should equal to `msstage`."))
        for ab in abundance
            length(chemicals) == length(ab) || throw(ArgumentError("The length of chemicals and each abundance are not equaled."))
        end
        id = sortperm(abundance[end]; rev = true)
        ab = hcat([ab[id] for ab in abundance]...)
        new{T, N}(chemicals[id], ab)
    end
    function Isobars(chemicals::Vector{T}, abundance::Matrix{N}) where {T<:ChemicalTransition, N}
        allequal(msstage, chemicals) || throw(ArgumentError("All chemicals should have the same `msstage`."))
        size(abundance, 1) == length(chemicals) || throw(ArgumentError("The column size of `abundance` should equal to the length of `chemicals`."))
        size(abundance, 2) == msstage(chemicals[begin]) || throw(ArgumentError("The row size of `abundance` should equal to `msstage`."))
        id = sortperm(abundance[:, end]; rev = true)
        new{T, N}(chemicals[id], abundance[id, :])
    end
end

Isobars(chemicals::AbstractVector, abundance::AbstractArray) = _Isobars(collect(chemicals), abundance)
_Isobars(chemicals::AbstractVector, abundance::AbstractArray) = Isobars(chemicals, collect(abundance))
_Isobars(chemicals::AbstractVector, abundance::VecOrMat) = Isobars(chemicals, abundance)

struct ElementsVector
    elements::Vector{String}
    numbers::Vector{Int}
end

"""
    Isotopomers{T<:AbstractChemical} <: AbstractChemical

Chemicals differed from isotopic replacement location.

# Fields 
* `parent::T`: shared chemical structure of isotopomers prior to isotopic replacement. 
* `isotopes::ElementsVector`: isotopes-number pairs of isotopic replacement.

# Constructors
    Isotopomers(parent::AbstractChemical, isotopes::ElementsVector)
    Isotopomers(parent::AbstractChemical, fullformula::String)
    Isotopomers(parent::AbstractChemical, fullelements::Dict)
    Isotopomers(parent::AbstractChemical, fullelements::Vector{Pair{String, Int}})

All minor isotopes are regarded as isotopic replacement in `fullformula` and `fullelements`.
"""
struct Isotopomers{T<:AbstractChemical} <: AbstractChemical
    parent::T 
    isotopes::ElementsVector
end

function Isotopomers(chemical::AbstractChemical, fullformula::String)
    Isotopomers(chemicalparent(chemical), dictionary_elements(chemicalelements(fullformula)))
end

function Isotopomers(chemical::AbstractChemical, fullelements::Vector{Pair{String, Int}})
    Isotopomers(chemicalparent(chemical), dictionary_elements(fullelements))
end

function Isotopomers(chemical::AbstractChemical, fullelements::Dict)
    parent = chemicalparent(chemical)
    dp = dictionary_elements(chemicalelements(parent))
    ev = ElementsVector(collect(keys(fullelements)), collect(values(fullelements)))
    del = Int[]
    for (i, (k, n)) in enumerate(ev)
        iselement(k) && (push!(del, i); continue)
        ev.numbers[i] = n - get(dp, k, 0) 
    end
    deleteat!(ev.elements, del)
    deleteat!(ev.numbers, del)
    Isotopomers(parent, ev)
end

"""
    Groupedisotopomers{T<:AbstractChemical, N} <: AbstractChemical

Isotopomerized chemicals grouped by isotopomer state.

# Fields 
* `parent::T`: shared chemical structure prior to isotopic replacement. 
* `state::Int`: isotopomer state.
* `isotope::String`: isotope for computing isotopomer state.
* `isotopes::Vector{ElementsVector}`: Isotopes-number pairs of isotopic replacements of each isotopomers.
* `abundance::Vector{N}`: abundance of each isotopomers.
"""
struct Groupedisotopomers{T<:AbstractChemical, N} <: AbstractChemical
    parent::T 
    state::Int
    isotope::String
    isotopes::Vector{ElementsVector}
    abundance::Vector{N}
    function Groupedisotopomers(parent::T, state::Int, isotope::String, isotopes::Vector{ElementsVector}, abundance::Vector{N}) where {T, N}
        id = sortperm(abundance)
        new{T, N}(parent, state, isotope, isotopes[id], abundance[id])
    end
end

"""
    ChemicalSeries(chemical::AbstractChemicalsSchema)
    ChemicalSeries(pair::Pair)
    ChemicalSeries(chemicals::AbstractVector)

Transform chemical into valid chemical structure. Multiple chemicals are converted into `ChemicalTransition`.
"""
ChemicalSeries(cc::AbstractChemicalsSchema) = cc
ChemicalSeries(cc::ChemicalTransition) = cc
ChemicalSeries(ct...) = ChemicalTransition(ct...) 
ChemicalSeries(v::AbstractVector) = length(v) < 2 ? ChemicalSeries(first(v)) : ChemicalTransition(v...)
ChemicalSeries(v::Pair) = ChemicalTransition(_ChemicalSeries(v)...)
_ChemicalSeries(v::Pair) = (_ChemicalSeries(first(v))..., _ChemicalSeries(last(v))...)
_ChemicalSeries(v::ChemicalTransition) = (chemicaltransition(v)..., )
_ChemicalSeries(v::AbstractChemicalsSchema) = (v, ) 

"""
    AbstractChemicalWrapper{T<:AbstractChemical} <: AbstractChemical 
    
Abstract type for all types wrapping a chemical of type `T` as field `chemical`. By default, all attributes come from `chemical`. 
"""
abstract type AbstractChemicalWrapper{T<:AbstractChemical} <: AbstractChemical end

"""
    Electron{T} <: AbstractChemicalWrapper{T}
    Electron(constructor, args...; kwargs...) 
    Electron(name = "Electron", elements = Pair{String, Int}[]; charge = -1, abbreviation = "e")

Electron.

# Fields 
* `chemical::T`.
"""
struct Electron{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Proton{T} <: AbstractChemicalWrapper{T}
    Proton(constructor, args...; kwargs...) 
    Proton(name = "Proton", elements = ["H" => 1]; charge = 1, abbreviation = "H")

Proton.

# Fields 
* `chemical::T`.
"""
struct Proton{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Water{T} <: AbstractChemicalWrapper{T}
    Water(constructor, args...; kwargs...) 
    Water(name = "Water", elements = ["H" => 2, "O" => 1]; charge = 0, abbreviation = "H2O")

Water.

# Fields 
* `chemical::T`.

# Attributes (default)
* `name`: `"Water"`.
* `chemicalelements`: `Pair{String, Int}["H" => 2, "O" => 1]`.
* `charge`: `0`.
* `abbreviation`: `"H2O"`.
"""
struct Water{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Ammonia{T} <: AbstractChemicalWrapper{T}
    Ammonia(constructor, args...; kwargs...) 
    Ammonia(name = "Ammonia", elements = ["N" => 1, "H" => 3]; charge = 0, abbreviation = "NH3")

Ammonia.

# Fields 
* `chemical::T`.
"""
struct Ammonia{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Ammonium{T} <: AbstractChemicalWrapper{T}
    Ammonium(constructor, args...; kwargs...) 
    Ammonium(name = "Ammonium", elements = ["N" => 1, "H" => 4]; charge = 1, abbreviation = "NH4")

Ammonium.

# Fields 
* `chemical::T`.
"""
struct Ammonium{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Sodium{T} <: AbstractChemicalWrapper{T}
    Sodium(constructor, args...; kwargs...) 
    Sodium(name = "Sodium", elements = ["Na" => 1]; charge = 1, abbreviation = "Na")

Sodium.

# Fields 
* `chemical::T`.
"""
struct Sodium{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Potassium{T} <: AbstractChemicalWrapper{T}
    Potassium(constructor, args...; kwargs...) 
    Potassium(name = "Potassium", elements = ["K" => 1]; charge = 1, abbreviation = "K")

Potassium.

# Fields 
* `chemical::T`.
"""
struct Potassium{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Lithium{T} <: AbstractChemicalWrapper{T}
    Lithium(constructor, args...; kwargs...) 
    Lithium(name = "Lithium", elements = ["Li" => 1]; charge = 1, abbreviation = "Li")

Lithium.

# Fields 
* `chemical::T`.
"""
struct Lithium{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Silver{T} <: AbstractChemicalWrapper{T}
    Silver(constructor, args...; kwargs...) 
    Silver(name = "Silver", elements = ["Ag" => 1]; charge = 1, abbreviation = "Ag")

Silver.

# Fields 
* `chemical::T`.
"""
struct Silver{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Acetate{T} <: AbstractChemicalWrapper{T}
    Acetate(constructor, args...; kwargs...) 
    Acetate(name = "Acetate", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1]; 
            charge = -1, abbreviation = "OAc")

Acetate.

# Fields 
* `chemical::T`.
"""
struct Acetate{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Formate{T} <: AbstractChemicalWrapper{T}
    Formate(constructor, args...; kwargs...) 
    Formate(name = "Formate", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1]; 
            charge = -1, abbreviation = "OFo")

Formate.

# Fields 
* `chemical::T`.
"""
struct Formate{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    AceticAcid{T} <: AbstractChemicalWrapper{T}
    AceticAcid(constructor, args...; kwargs...) 
    AceticAcid(name = "AceticAcid", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1]; 
            charge = 0, abbreviation = "HOAc")

Acetic acid.

# Fields 
* `chemical::T`.
"""
struct AceticAcid{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    FormicAcid{T} <: AbstractChemicalWrapper{T}
    FormicAcid(constructor, args...; kwargs...) 
    FormicAcid(name = "FormicAcid", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]; 
            charge = 0, abbreviation = "HOFo")

Formic acid.

# Fields 
* `chemical::T`.
"""
struct FormicAcid{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    MethylAcetate{T} <: AbstractChemicalWrapper{T}
    MethylAcetate(constructor, args...; kwargs...) 
    MethylAcetate(name = "MethylAcetate", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3]; 
            charge = 0, abbreviation = "MeOAc")

Methyl acetate.

# Fields 
* `chemical::T`.
"""
struct MethylAcetate{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    MethylFormate{T} <: AbstractChemicalWrapper{T}
    MethylFormate(constructor, args...; kwargs...) 
    MethylFormate(name = "MethylFormate", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3]; 
            charge = 0, abbreviation = "MeOFo")

Methyl formate.

# Fields 
* `chemical::T`.
"""
struct MethylFormate{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Fluoride{T} <: AbstractChemicalWrapper{T}
    Fluoride(constructor, args...; kwargs...) 
    Fluoride(name = "Fluoride", elements = ["F" => 1]; charge = -1, abbreviation = "F")

Fluoride.

# Fields 
* `chemical::T`.
"""
struct Fluoride{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Chloride{T} <: AbstractChemicalWrapper{T}
    Chloride(constructor, args...; kwargs...) 
    Chloride(name = "Chloride", elements = ["Cl" => 1]; charge = -1, abbreviation = "Cl")

Chloride.

# Fields 
* `chemical::T`.
"""
struct Chloride{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Methenium{T} <: AbstractChemicalWrapper{T}
    Methenium(constructor, args...; kwargs...) 
    Methenium(name = "Methenium", elements = ["C" => 1, "H" => 3]; charge = -1, abbreviation = "Me")

Methenium.

# Fields 
* `chemical::T`.

# Attributes (default)
* `name`: `"Methenium"`.
* `chemicalelements`: `Pair{String, Int}["C" => 1, "H" => 3]`.
* `charge`: `1`.
* `abbreviation`: `"Me"`.
"""
struct Methenium{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 

for fn in ["Electron", "Proton", "Water", "Ammonia", "Ammonium", "Sodium", "Potassium", "Lithium", "Silver", "Acetate", "Formate", "AceticAcid", "FormicAcid", "MethylAcetate", "MethylFormate", "Fluoride", "Chloride", "Methenium"]
    expr_l = :(fn(constructor::Type{<:AbstractChemical}, args...; kwargs...)) 
    expr_r = :(fn(constructor(args...; kwargs...)))
    expr_l.args[1] = Symbol(fn)
    expr_r.args[1] = Symbol(fn)
    eval(Expr(:(=), expr_l, expr_r))
end

Electron(name::AbstractString = "Electron", elements = Pair{String, Int}[]; charge = -1, abbreviation = "e", kwargs...) = Electron(Chemical(name, elements; charge, abbreviation, kwargs...))
Proton(name::AbstractString = "Proton", elements = ["H" => 1]; charge = 1, abbreviation = "H", kwargs...) = Proton(Chemical(name, elements; charge, abbreviation, kwargs...))
Water(name::AbstractString = "Water", elements = ["H" => 2, "O" => 1]; charge = 0, abbreviation = "H2O", kwargs...) = Water(Chemical(name, elements; charge, abbreviation, kwargs...))
Ammonia(name::AbstractString = "Ammonia", elements = ["N" => 1, "H" => 3]; charge = 0, abbreviation = "NH3", kwargs...) = Ammonia(Chemical(name, elements; charge, abbreviation, kwargs...))
Ammonium(name::AbstractString = "Ammonium", elements = ["N" => 1, "H" => 4]; charge = 1, abbreviation = "NH4", kwargs...) = Ammonium(Chemical(name, elements; charge, abbreviation, kwargs...))
Sodium(name::AbstractString = "Sodium", elements = ["Na" => 1]; charge = 1, abbreviation = "Na", kwargs...) = Sodium(Chemical(name, elements; charge, abbreviation, kwargs...))
Potassium(name::AbstractString = "Potassium", elements = ["K" => 1]; charge = 1, abbreviation = "K", kwargs...) = Potassium(Chemical(name, elements; charge, abbreviation, kwargs...))
Lithium(name::AbstractString = "Lithium", elements = ["Li" => 1]; charge = 1, abbreviation = "Li", kwargs...) = Lithium(Chemical(name, elements; charge, abbreviation, kwargs...))
Silver(name::AbstractString = "Silver", elements = ["Ag" => 1]; charge = 1, abbreviation = "Ag", kwargs...) = Silver(Chemical(name, elements; charge, abbreviation, kwargs...))
Acetate(name::AbstractString = "Acetate", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1]; charge = -1, abbreviation = "OAc", kwargs...) = Acetate(Chemical(name, elements; charge, abbreviation, kwargs...))
Formate(name::AbstractString = "Formate", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1]; charge = -1, abbreviation = "OFo", kwargs...) = Formate(Chemical(name, elements; charge, abbreviation, kwargs...))
AceticAcid(name::AbstractString = "Acetic Acid", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1]; charge = 0, abbreviation = "HOAc", kwargs...) = AceticAcid(Chemical(name, elements; charge, abbreviation, kwargs...))
FormicAcid(name::AbstractString = "Formic Acid", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]; charge = -1, abbreviation = "HOFo", kwargs...) = FormicAcid(Chemical(name, elements; charge, abbreviation, kwargs...))
MethylAcetate(name::AbstractString = "Methyl Acetate", elements = ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3]; charge = 0, abbreviation = "MeOAc", kwargs...) = MethylAcetate(Chemical(name, elements; charge, abbreviation, kwargs...))
MethylFormate(name::AbstractString = "Methyl Formate", elements = ["H" => 1, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3]; charge = 0, abbreviation = "MeOFo", kwargs...) = MethylFormate(Chemical(name, elements; charge, abbreviation, kwargs...))
Fluoride(name::AbstractString = "Fluoride", elements = ["F" => 1]; charge = -1, abbreviation = "F", kwargs...) = Fluoride(Chemical(name, elements; charge, abbreviation, kwargs...))
Chloride(name::AbstractString = "Chloride", elements = ["Cl" => 1]; charge = -1, abbreviation = "Cl", kwargs...) = Chloride(Chemical(name, elements; charge, abbreviation, kwargs...))
Methenium(name::AbstractString = "Methenium", elements = ["C" => 1, "H" => 3]; charge = -1, abbreviation = "Me", kwargs...) = Methenium(Chemical(name, elements; charge, abbreviation, kwargs...))

"""
    const GenericChemical = Union{Chemical, FormulaChemical, <:AbstractChemicalWrapper{Chemical}, <:AbstractChemicalWrapper{FormulaChemical}}

Generic chemical types.
"""
const GenericChemical = Union{Chemical, FormulaChemical, <:AbstractChemicalWrapper{Chemical}, <:AbstractChemicalWrapper{FormulaChemical}}

"""
    AbstractAdductIon{S, T} <: AbstractChemical

Abstract type for adduct ions with core chemical type `S` and adduct type `T`.

# Special attributes
* [`ncore`](@ref) `-> Int`: number of core chemical "M" in adduct ion representation "[M+X]n+". 
* [`ioncore`](@ref) `-> S`: the core chemical undergoing ionization. 
* [`ionadduct`](@ref) `-> T`: the adduct formed during ionization. 
"""
abstract type AbstractAdductIon{S, T} <: AbstractChemical end

"""
    AdductIon{S<:AbstractChemical, T<:AbstractScheme} <: AbstractAdductIon{S, T}

Adduct ions formed in mass spectrometry.

# Fields
* `core`: the core chemical undergoing ionization. 
* `adduct`: the adduct formed during ionization.
* `ncore`: number of core chemical undergoing ionization. 

# Constructors
    AdductIon(core::AbstractChemical, adduct::AbstractScheme, ncore = 1)
    AdductIon(core::AbstractChemical, adduct_string::AbstractString)

`adduct_string` is parsed by `parse_adduct(adduct_string; args = true)`.
"""
struct AdductIon{S<:AbstractChemical, T<:AbstractScheme} <: AbstractAdductIon{S, T}
    core::S
    adduct::T
    ncore::Int
    function AdductIon(core::S, adduct::AbstractScheme, ncore::Int) where {S<:AbstractChemical}
        adduct = completescheme(core, adduct)
        new{S, typeof(adduct)}(core, adduct, ncore)
    end
end

AdductIon(cc::AbstractChemical, a::AbstractScheme) = AdductIon(cc, a, 1)
AdductIon(cc::AbstractChemical, a::AbstractString, n = 1) = AdductIon(cc, parse_adduct(a; args = true)...)