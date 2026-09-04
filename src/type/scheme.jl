"""
    AbstractElementalScheme <: AbstractScheme end

Abstract scheme contaning exact elements.
"""
abstract type AbstractElementalScheme <: AbstractScheme end
"""
    AbstractStructuralScheme <: AbstractScheme end

Abstract scheme contaning only structures.  
"""
abstract type AbstractStructuralScheme <: AbstractScheme end
"""
    AbstractCompleteScheme{T, S} <: AbstractScheme end

Abstract scheme contaning both elements and structures. 
"""
abstract type AbstractCompleteScheme{T, S} <: AbstractScheme end
"""
    StructuralChemicalScheme <: AbstractStructuralScheme end

Abstract scheme contaning structures generating chemical entity. 
"""
abstract type StructuralChemicalScheme <: AbstractStructuralScheme end

struct RandomProductScheme <: AbstractStructuralScheme end

"""
    StructuralElementalScheme{T, S} <: AbstractCompleteScheme{T, S} end

Default [`AbstractCompleteScheme`](@ref).

# Fields
* `structuralscheme::T`.
* `elementalscheme::S`.
"""
struct StructuralElementalScheme{T, S} <: AbstractCompleteScheme{T, S}
    structuralscheme::T 
    elementalscheme::S
end

"""

    ElementalScheme{Bool, T<:AbstractChemical} <: AbstractElementalScheme
    ElementalScheme(gain::Bool, chemical::T) = ElementalScheme{gain, T}(chemical)
    
Single scheme involving a chemical. The elements are fixed, and can be replaced by minor isotopes. 
* `ElementalScheme{false}`: chemical loss from a precursor. This product is not detected in MS; the other part of precursor is detected instead.
* `ElementalScheme{true}`: chemical gain to a precursor. This product is not detected in MS; the merged chemical is detected instead.

# Fields 
* `chemical::T`: chemical involved in scheme.

Single isotopomer can be set by using [`Isotopomers`](@ref) as field `chemical`. 
Elemental scheme can be redirected to the corresponding isotopic labeled scheme in [`AdductIon`](@ref) by dispatching on core chemical and existing schema, 
or looking up the `property` for generic [`Chemical`](@ref).
"""
struct ElementalScheme{Bool, T<:AbstractChemical} <: AbstractElementalScheme
    chemical::T
    function ElementalScheme(gain::Bool, x::T) where {T<:AbstractChemical}
        new{gain, T}(x)
    end
end

"""
    ChemicalGain(chemical)

Chemical gain of `chemical`, i.e. `ElementalScheme(true, chemical)`. See [`ElementalScheme`](@ref).
"""
ChemicalGain(x) = ElementalScheme(true, x)

"""
    ChemicalLoss(chemical)

Chemical loss of `chemical`, i.e. `ElementalScheme(false, chemical)`. See [`ElementalScheme`](@ref).
"""
ChemicalLoss(x) = ElementalScheme(false, x)

"""
    ChemicalSchema{T<:AbstractScheme} <: AbstractScheme

Mutiple chemical schema.

# Fields
* `schema::Vector{Pair{T, Int}}`: scheme => number vector. The number v is the times of the scheme in the chemical schema. 
"""
struct ChemicalSchema{T<:AbstractScheme} <: AbstractScheme
    schema::Vector{T}
    number::Vector{Int}
end 

"""
    IsotopomerizedSchema{T<:ChemicalSchema} <: AbstractScheme

Mutiple chemical schema with delocalized isotopic replacements.

# Fields
* `schema::T`: parent scheme.
* `isotopes::ElementsVector`: delocalized isotopic replacements.
"""
struct IsotopomerizedSchema{T<:ChemicalSchema} <: AbstractScheme 
    parent::T
    isotopes::ElementsVector
end

function IsotopomerizedSchema(chemical::AbstractScheme, fullformula::String)
    IsotopomerizedSchema(chemicalparent(chemical), dictionary_elements(chemicalelements(fullformula)))
end

function IsotopomerizedSchema(chemical::AbstractScheme, fullelements::Vector{Pair{String, Int}})
    IsotopomerizedSchema(chemicalparent(chemical), dictionary_elements(fullelements))
end

function IsotopomerizedSchema(chemical::AbstractScheme, fullelements::Dict)
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
    IsotopomerizedSchema(parent, ev)
end


"""
    Groupedisotopomerizedschema{T<:AbstractScheme, N} <: AbstractScheme

Isotopomerized schema grouped by isotopomer state.

# Fields 
* `parent::T`: shared chemical scheme prior to isotopic replacement. 
* `state::Int`: isotopomer state.
* `isotope::String`: isotope for computing isotopomer state.
* `isotopes::Vector{ElementsVector}`: Isotopes-number pairs of isotopic replacements of each isotopomers.
* `abundance::Vector{N}`: abundance of each isotopomers.
"""
struct Groupedisotopomerizedschema{T<:AbstractScheme, N} <: AbstractScheme
    parent::T 
    state::Int
    isotope::String
    isotopes::Vector{ElementsVector}
    abundance::Vector{N}
    function Groupedisotopomerizedschema(parent::T, state::Int, isotope::String, isotopes::Vector{ElementsVector}, abundance::Vector{N}) where {T, N}
        id = sortperm(abundance)
        new{T, N}(parent, state, isotope, isotopes[id], abundance[id])
    end
end

schemetype(::ChemicalSchema{T}) where T = T 
schemetype(::T) where T = T 

function ChemicalSchema(scheme::T, schema...) where {T<:AbstractScheme} 
    C = promote_type(T, schemetype.(schema)...)
    cs = C[scheme]
    cn = Int[1]
    for s in schema
        push_scheme!(cs, cn, s)
    end
    ChemicalSchema(cs, cn)
end

function ChemicalSchema(schema::AbstractVector{T}) where {T<:AbstractScheme}
    cs = T[first(schema)]
    cn = Int[1]
    length(schema) < 2 && return ChemicalSchema(cs, cn)
    for s in @view schema[2:end]
        push_scheme!(cs, cn, s)
    end
    ChemicalSchema(cs, cn)
end

function ChemicalSchema(scheme::ChemicalSchema{T}, schema...) where {T<:AbstractScheme}
    C = promote_type(T, schemetype.(schema)...)
    if C == T
        cs = copy(scheme.schema)
    else
        cs = convert(Vector{C}, copy(scheme.schema))
    end
    cn = copy(scheme.number)
    for s in schema
        push_scheme!(cs, cn, s)
    end
    ChemicalSchema(cs, cn)
end

function push_scheme!(cs::Vector, cn::Vector, scheme::ChemicalSchema)
    for (k, v) in zip(scheme.schema, scheme.number)
        i = findfirst(==(k), cs)
        if i !== nothing
            cn[i] += v
        else
            push!(cs, k)
            push!(cn, v)
        end
    end
    cs
end

function push_scheme!(cs::Vector, cn::Vector, scheme::AbstractScheme)
    i = findfirst(==(scheme), cs)
    if i !== nothing
        cn[i] += 1
    else
        push!(cs, scheme)
        push!(cn, 1)
    end
    cs
end

"""
    const CompleteSchema = Union{<:AbstractCompleteScheme, <:ChemicalSchema{<:AbstractCompleteScheme}, <:IsotopomerizedSchema{<:ChemicalSchema{<:AbstractCompleteScheme}}, <:Groupedisotopomerizedschema{<:ChemicalSchema{<:AbstractCompleteScheme}}}

Complete scheme (scheme containing both [`structuralscheme`](@ref) and [`elementalscheme`](@ref)).
"""
const CompleteSchema = Union{<:AbstractCompleteScheme, <:ChemicalSchema{<:AbstractCompleteScheme}, <:IsotopomerizedSchema{<:ChemicalSchema{<:AbstractCompleteScheme}}, <:Groupedisotopomerizedschema{<:ChemicalSchema{<:AbstractCompleteScheme}}}

"""
    const StructuralSchema = Union{<:AbstractStructuralScheme, <:ChemicalSchema{<:AbstractStructuralScheme}, <:IsotopomerizedSchema{<:ChemicalSchema{<:AbstractStructuralScheme}}, <:Groupedisotopomerizedschema{<:ChemicalSchema{<:AbstractStructuralScheme}}}

Stuctural schema.
"""
const StructuralSchema = Union{<:AbstractStructuralScheme, <:ChemicalSchema{<:AbstractStructuralScheme}, <:IsotopomerizedSchema{<:ChemicalSchema{<:AbstractStructuralScheme}}, <:Groupedisotopomerizedschema{<:ChemicalSchema{<:AbstractStructuralScheme}}}

"""
    const ElementalSchema = Union{<:AbstractElementalScheme, <:ChemicalSchema{<:AbstractElementalScheme}, <:IsotopomerizedSchema{<:ChemicalSchema{<:AbstractElementalScheme}}, <:Groupedisotopomerizedschema{<:ChemicalSchema{<:AbstractElementalScheme}}}

Elemental schema.
"""
const ElementalSchema = Union{<:AbstractElementalScheme, <:ChemicalSchema{<:AbstractElementalScheme}, <:IsotopomerizedSchema{<:ChemicalSchema{<:AbstractElementalScheme}}, <:Groupedisotopomerizedschema{<:ChemicalSchema{<:AbstractElementalScheme}}}

"""
    const CompleteSchemeChemical = AbstractCompleteScheme{T, <:AbstractChemical}

Complete scheme genrating chemical.
"""
const CompleteSchemeChemical = AbstractCompleteScheme{T, <:AbstractChemical} where T