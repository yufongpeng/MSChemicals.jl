"""
    Chemical <: AbstractChemical

Unstructured chemical type with its name, elements (formula), and additional properties.

# Fields 
* `name::String`: a unique chemical name.
* `elements::Vector{Pair{String, Int}}`: chemical elements.
* `property::Vector{Pair{Symbol, Any}}`: additional properties; the pairs repressent names and values.

# Constructors
* `Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}, property::Vector{Pair{Symbol, Any}})`
* `Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}; kwargs_as_property_pairs...)`
* `Chemical(name::AbstractString, formula::AbstractString; kwargs_as_property_pairs...)`
"""
struct Chemical <: AbstractChemical
    name::String
    elements::Vector{Pair{String, Int}}
    property::Vector{Pair{Symbol, Any}}
end

Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}; kwargs...) = Chemical(name, elements, collect(kwargs))
Chemical(name::AbstractString, formula::AbstractString; kwargs...) = Chemical(name, chemicalelements(formula), collect(kwargs))

"""
    FormulaChemical <: AbstractChemical

Unstructured chemical type with elements (formula), and additional properties. Chemical name will be formula.

# Fields 
* `elements::Vector{Pair{String, Int}}`: chemical elements.
* `property::Vector{Pair{Symbol, Any}}`: additional properties; the pairs repressent names and values.

# Constructors
* `FormulaChemical(elements::Vector{Pair{String, Int}}, property::Vector{Pair{Symbol, Any}})`
* `FormulaChemical(elements::Vector{Pair{String, Int}}; kwargs_as_property_pairs...)`
* `FormulaChemical(formula::AbstractString; kwargs_as_property_pairs...)`
"""
struct FormulaChemical <: AbstractChemical
    elements::Vector{Pair{String, Int}}
    property::Vector{Pair{Symbol, Any}}
end

FormulaChemical(elements::Vector{Pair{String, Int}}; kwargs...) = FormulaChemical(elements, collect(kwargs))
FormulaChemical(formula::AbstractString; kwargs...) = FormulaChemical(chemicalelements(formula), collect(kwargs))

"""
    ChemicalTransition{T<:AbstractChemicalsSchema} <: AbstractChemical

Chemical transition in MSⁿ. Products can be a `ChemicalLoss` or `ChemicalGain`.

# Fields 
* `transition::Vector{T}`

# Constructors
* `ChemicalTransition(transition::Vector)`
* `ChemicalTransition(precursor, products...)`: push `products` into `precursor`.
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
* `Isobars(chemicals::Vector, abundance::Vector)`
* `Isobars(chemicals::Vector{<: ChemicalTransition}, abundance::Vector)`
* `Isobars(chemicals::Vector{<: ChemicalTransition}, abundance::Matrix)`
"""
struct Isobars{T<:AbstractChemical, N} <: AbstractChemical
    chemicals::Vector{T}
    abundance::VecOrMat{N}
    function Isobars(chemicals::Vector{T}, abundance::Vector{N}) where {T, N}
        id = sortperm(abundance; rev = true)
        new{T, N}(chemicals[id], abundance[id])
    end
    function Isobars(chemicals::Vector{T}, abundance::Vector{N}) where {T <: ChemicalTransition, N}
        allequal(msstage, chemicals) || throw(ArgumentError("All chemicals should have the same `msstage`."))
        id = sortperm(abundance; rev = true)
        ab = hcat([abundance[id] for _ in eachindex(chemicals)]...)
        new{T, N}(chemicals[id], ab)
    end
    function Isobars(chemicals::Vector{T}, abundance::Matrix{N}) where {T <: ChemicalTransition, N}
        allequal(msstage, chemicals) || throw(ArgumentError("All chemicals should have the same `msstage`."))
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
* `Isotopomers(parent::AbstractChemical, isotopes::ElementsVector)`
* `Isotopomers(parent::AbstractChemical, fullformula::String)`
* `Isotopomers(parent::AbstractChemical, fullelements::Dict)`
* `Isotopomers(parent::AbstractChemical, fullelements::Vector{Pair{String, Int}})`
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
    AbstractChemicalWrapper{T <: AbstractChemical} <: AbstractChemical 
    
Abstract type for all types wrapping a chemical of type `T` as field `chemical`. By default, all attributes come from `chemical`. 
"""
abstract type AbstractChemicalWrapper{T <: AbstractChemical} <: AbstractChemical end

"""
    Electron{T} <: AbstractChemicalWrapper{T}
    Electron()

Electron

# Attributes
* `name`: `"Electron"`
* `chemicalelements`: `Pair{String, Int}[]`
* `charge`: `-1`
* `abbreviation`: `"e"`
"""
struct Electron{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Proton{T} <: AbstractChemicalWrapper{T}
    Proton()

Proton

# Attributes
* `name`: `"Proton"`
* `chemicalelements`: `Pair{String, Int}["H" => 1]`
* `charge`: `1`
* `abbreviation`: `"H"`
"""
struct Proton{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Water{T} <: AbstractChemicalWrapper{T}
    Water()

Water

# Attributes
* `name`: `"Water"`
* `chemicalelements`: `Pair{String, Int}["H" => 2, "O" => 1]`
* `charge`: `0`
* `abbreviation`: `"H2O"`
"""
struct Water{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Ammonia{T} <: AbstractChemicalWrapper{T}
    Ammonia()

Ammonia

# Attributes
* `name`: `"Ammonia"`
* `chemicalelements`: `Pair{String, Int}["N" => 1, "H" => 3]`
* `charge`: `0`
* `abbreviation`: `"NH3"`
"""
struct Ammonia{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Ammonium{T} <: AbstractChemicalWrapper{T}
    Ammonium()

Ammonium

# Attributes
* `name`: `"Ammonium"`
* `chemicalelements`: `Pair{String, Int}["N" => 1, "H" => 4]`
* `charge`: `1`
* `abbreviation`: `"NH4"`
"""
struct Ammonium{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Sodium{T} <: AbstractChemicalWrapper{T}
    Sodium()

Sodium

# Attributes
* `name`: `"Sodium"`
* `chemicalelements`: `Pair{String, Int}["Na" => 1]`
* `charge`: `1`
* `abbreviation`: `"Na"`
"""
struct Sodium{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Potassium{T} <: AbstractChemicalWrapper{T}
    Potassium()

Potassium

# Attributes
* `name`: `"Potassium"`
* `chemicalelements`: `Pair{String, Int}["K" => 1]`
* `charge`: `1`
* `abbreviation`: `"K"`
"""
struct Potassium{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Lithium{T} <: AbstractChemicalWrapper{T}
    Lithium()

Lithium

# Attributes
* `name`: `"Lithium"`
* `chemicalelements`: `Pair{String, Int}["Li" => 1]`
* `charge`: `1`
* `abbreviation`: `"Li"`
"""
struct Lithium{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Silver{T} <: AbstractChemicalWrapper{T}
    Silver()

Silver

# Attributes
* `name`: `"Silver"`
* `chemicalelements`: `Pair{String, Int}["Ag" => 1]`
* `charge`: `1`
* `abbreviation`: `"Ag"`
"""
struct Silver{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Acetate{T} <: AbstractChemicalWrapper{T}
    Acetate()

Acetate

# Attributes
* `name`: `"Acetate"`
* `chemicalelements`: `Pair{String, Int}["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1]`
* `charge`: `-1`
* `abbreviation`: `"OAc"`
"""
struct Acetate{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Formate{T} <: AbstractChemicalWrapper{T}
    Formate()

Formate

# Attributes
* `name`: `"Formate"`
* `chemicalelements`: `Pair{String, Int}["H" => 1, "C" => 1, "O" => 1, "O" => 1]`
* `charge`: `-1`
* `abbreviation`: `"OFo"`
"""
struct Formate{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    AceticAcid{T} <: AbstractChemicalWrapper{T}
    AceticAcid()

Acetic acid

# Attributes
* `name`: `"AceticAcid"`
* `chemicalelements`: `Pair{String, Int}["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1]`
* `charge`: `0`
* `abbreviation`: `"HOAc"`
"""
struct AceticAcid{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    FormicAcid{T} <: AbstractChemicalWrapper{T}
    FormicAcid()

Formic acid

# Attributes
* `name`: `"FormicAcid"`
* `chemicalelements`: `Pair{String, Int}["H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]`
* `charge`: `0`
* `abbreviation`: `"HOFo"`
"""
struct FormicAcid{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    MethylAcetate{T} <: AbstractChemicalWrapper{T}
    MethylAcetate()

Methyl acetate

# Attributes
* `name`: `"MethylAcetate"`
* `chemicalelements`: `Pair{String, Int}["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3]`
* `charge`: `0`
* `abbreviation`: `"MeOAc"`
"""
struct MethylAcetate{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    MethylFormate{T} <: AbstractChemicalWrapper{T}
    MethylFormate()

Methyl formate

# Attributes
* `name`: `"MethylFormate"`
* `chemicalelements`: `Pair{String, Int}["H" => 1, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3]`
* `charge`: `0`
* `abbreviation`: `"MeOFo"`
"""
struct MethylFormate{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Fluoride{T} <: AbstractChemicalWrapper{T}
    Fluoride()

Fluoride

# Attributes
* `name`: `"Fluoride"`
* `chemicalelements`: `Pair{String, Int}["F" => 1]`
* `charge`: `-1`
* `abbreviation`: `"F"`
"""
struct Fluoride{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Chloride{T} <: AbstractChemicalWrapper{T}
    Chloride()

Chloride

# Attributes
* `name`: `"Chloride"`
* `chemicalelements`: `Pair{String, Int}["Cl" => 1]`
* `charge`: `-1`
* `abbreviation`: `"Cl"`
"""
struct Chloride{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 
"""
    Methenium{T} <: AbstractChemicalWrapper{T}
    Methenium()

Methenium

# Attributes
* `name`: `"Methenium"`
* `chemicalelements`: `Pair{String, Int}["C" => 1, "H" => 3]`
* `charge`: `1`
* `abbreviation`: `"Me"`
"""
struct Methenium{T} <: AbstractChemicalWrapper{T}
    chemical::T
end 

Electron() = Electron(Chemical("Electron", Pair{String, Int}[]; charge = -1, abbreviation = "e"))
Proton() = Proton(Chemical("Proton", ["H" => 1]; charge = 1, abbreviation = "H"))
Water() = Water(Chemical("Water", ["H" => 2, "O" => 1]; charge = 0, abbreviation = "H2O"))
Ammonia() = Ammonia(Chemical("Ammonia", ["N" => 1, "H" => 3]; charge = 0, abbreviation = "NH3"))
Ammonium() = Ammonium(Chemical("Ammonium", ["N" => 1, "H" => 4]; charge = 1, abbreviation = "NH4"))
Sodium() = Sodium(Chemical("Sodium", ["Na" => 1]; charge = 1, abbreviation = "Na"))
Potassium() = Potassium(Chemical("Potassium", ["K" => 1]; charge = 1, abbreviation = "K"))
Lithium() = Lithium(Chemical("Lithium", ["Li" => 1]; charge = 1, abbreviation = "Li"))
Silver() = Silver(Chemical("Silver", ["Ag" => 1]; charge = 1, abbreviation = "Ag"))
Acetate() = Acetate(Chemical("Acetate", ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1]; charge = -1, abbreviation = "OAc"))
Formate() = Formate(Chemical("Formate", ["H" => 1, "C" => 1, "O" => 1, "O" => 1]; charge = -1, abbreviation = "OFo"))
AceticAcid() = AceticAcid(Chemical("Acetic Acid", ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "H" => 1]; charge = 0, abbreviation = "HOAc"))
FormicAcid() = FormicAcid(Chemical("Formic Acid", ["H" => 1, "C" => 1, "O" => 1, "O" => 1, "H" => 1]; charge = -1, abbreviation = "HOFo"))
MethylAcetate() = MethylAcetate(Chemical("Methyl Acetate", ["C" => 1, "H" => 3, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3]; charge = 0, abbreviation = "MeOAc"))
MethylFormate() = MethylFormate(Chemical("Methyl Formate", ["H" => 1, "C" => 1, "O" => 1, "O" => 1, "C" => 1, "H" => 3]; charge = 0, abbreviation = "MeOFo"))
Fluoride() = Fluoride(Chemical("Fluoride", ["F" => 1]; charge = -1, abbreviation = "F"))
Chloride() = Chloride(Chemical("Chloride", ["Cl" => 1]; charge = -1, abbreviation = "Cl"))
Methenium() = Methenium(Chemical("Methenium", ["C" => 1, "H" => 3]; charge = -1, abbreviation = "Me"))

"""
    const GenericChemical = Union{Chemical, FormulaChemical, <: AbstractChemicalWrapper{Chemical}, <: AbstractChemicalWrapper{FormulaChemical}}

Generic chemical types.
"""
const GenericChemical = Union{Chemical, FormulaChemical, <: AbstractChemicalWrapper{Chemical}, <: AbstractChemicalWrapper{FormulaChemical}}

"""
    AbstractAdductIon{S, T} <: AbstractChemical

Abstract type for adduct ions with core chemical type `S` and adduct type `T`.

# Special attributes
* `ncore -> Int`: number of core chemical "M" in adduct ion representation "[M+X]n+". 
* `ioncore -> S`: the core chemical undergoing ionization. 
* `ionadduct -> T`: the adduct formed during ionization. 
* `adductelements -> Vector{Pair{String, Int}}`: the elements changed with adduct of `adduct_ion`.
* `adductisotopes -> Vector{Pair{String, Int}}`: the elements changed when the core chemical has isotopic labeling that is lost in adduct formation. The returned vector is element-number pairs.
"""
abstract type AbstractAdductIon{S, T} <: AbstractChemical end

"""
    AdductIon{S <: AbstractChemical, T <: AbstractScheme} <: AbstractAdductIon{S, T}

Adduct ions forming in mass spectrometry.

# Fields
* `core`: the core chemical undergoing ionization. 
* `adduct`: the adduct formed during ionization.
* `ncore`: number of core chemical undergoing ionization. 

# Constructors
* `AdductIon(core::AbstractChemical, adduct::AbstractScheme, n = 1)` -> `AdductIon(core, adduct, n)`
* `AdductIon(core::AbstractChemical, adduct::AbstractString)` -> `AdductIon(core, parse_adduct(adduct)...)`
"""
struct AdductIon{S <: AbstractChemical, T <: AbstractScheme} <: AbstractAdductIon{S, T}
    core::S
    adduct::T
    ncore::Int
    function AdductIon(core::S, adduct::AbstractScheme, ncore::Int) where {S<:AbstractChemical}
        adduct = completescheme(core, adduct)
        new{S, typeof(adduct)}(core, adduct, ncore)
    end
end

AdductIon(cc::AbstractChemical, a::AbstractScheme) = AdductIon(cc, a, 1)
AdductIon(cc::AbstractChemical, a::AbstractString) = AdductIon(cc, parse_adduct(a; args = true)...)