"""
    getchemicalproperty(chemical::AbstractChemicalsSchema, property::Symbol, default::Nothing, type = Any) -> Any
    getchemicalproperty(chemical::AbstractChemicalsSchema, property::Symbol, default::T, type = T) -> T

Get `property` from `chemical`. 

This function defaults to finds the `property`, and returns `default` if it is not available. The return value is asserted to be type `type`.

For type [`Chemical`](@ref) and properties other than `name`, `elements` and `formula`, it iterates through `chemical.property`. 
If no matched property name is found, it returns `default`.

For type [`AbstractAdductIon`](@ref), it searches for properties of itself and then the properties of core chemical without specialized methods.
"""
getchemicalproperty(chemical::AbstractChemicalsSchema, property::Symbol, default::T) where T = getchemicalproperty(chemical, property, default, T)
getchemicalproperty(chemical::AbstractChemicalsSchema, property::Symbol, default::Nothing) = getchemicalproperty(chemical, property, default, Any)
getchemicalproperty(chemical::AbstractChemicalsSchema, property::Symbol, default, ::Type{T}) where T = hasproperty(chemical, property) ? getproperty(chemical, property)::T : default::T 
getchemicalproperty(adduct_ion::AbstractAdductIon, property::Symbol, default, type::Type{T})  where T = hasproperty(adduct_ion, property) ? getproperty(adduct_ion, property)::T : getchemicalproperty(ioncore(adduct_ion), property, default, type)
getchemicalproperty(chemical::AbstractChemicalWrapper, property::Symbol, default, ::Type{T}) where T = hasproperty(chemical, property) ? getproperty(chemical, property)::T : hasproperty(chemical.chemical, property) ? getproperty(chemical.chemical, property)::T : default::T 
function getchemicalproperty(chemical::Union{Chemical, FormulaChemical}, property::Symbol, default, ::Type{T}) where T
    hasproperty(chemical, property) && return getproperty(chemical, property)::T
    for (p, v) in chemical.property
        p == property && return v::T
    end
    default::T
end

"""
    chemicalname(chemical::AbstractChemicalsSchema; verbose = true, n = 1, loss = false, bracket = true, kwargs...) -> String

The name of `chemical`. It should be unique for each chemical object. 

# Generic Methods
* `AbstractChemicalsSchema`: property search.
* `AbstractAdductIon`: name composed of name of core chemical and adduct.
* `AbstractChemicalWrapper`: name of field `chemical`.

# Specific Methods 
* Species/Transition Level.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref).
    * property: `:name`.
    * default: [`defaultname`](@ref).

# Keyword arguments
* `verbose::Bool` determines whether includes all names or not for chemical species. If `verbose` is false, only the first (most abundant) chemical is included.
* `n::Int` determines number of chemical.
* `loss::Bool` determines whether the chemical is part of chemical loss, and sign flips are propagated into this chemical.
* `bracket::Bool` determines whether wrapping the output string with "[...]" and charge.
"""
function chemicalname(chemical::AbstractChemicalsSchema; n = 1, kwargs...) 
    nm = getchemicalproperty(chemical, :name, defaultname(chemical))
    n > 1 ? string(n, nm) : nm
end
chemicalname(chemical::AbstractChemicalWrapper; n = 1, kwargs...) = chemicalname(chemical.chemical; n, kwargs...)
function chemicalname(adduct_ion::AbstractAdductIon; loss = false, bracket = true, n = 1, kwargs...) 
    r = chemicalname(ioncore(adduct_ion); loss = bracket ? false : loss, bracket = false, n, kwargs...)
    a = chemicalabbr(ionadduct(adduct_ion); loss = bracket ? false : loss, bracket = false, n, kwargs...)
    n = ncore(adduct_ion)
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    if n > 1 
        r = string(n, r)
    end
    bracket ? string("[", r, a, "]", charge_repr(charge(adduct_ion))) : string(r, a)
end

"""
    chemicalformula(chemical::AbstractChemical; delim = "", unique = true, ischemical = true, loss = false, kwargs...) -> String
    chemicalformula(chemical::AbstractScheme; delim = "", unique = true, ischemical = false, loss = false, kwargs...) -> String

The formula of chemical entity of `chemical`.

# Generic Methods
* `AbstractChemical`: property search.
* `AbstractAdductIon`: formuala combining core chemical and adduct.
* `AbstractChemicalWrapper`: formula of field `chemical`.
* `AbstractScheme`: property search; started with `"+"` for chemical gain and `"-"` for chemical loss.
* `AbstractStructuralScheme`: throw error.

# Specific Methods 
* Entity Level.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref) -> `reverse_formula`.
    * property: `:formula`.
    * default: `""`.
2. Function: [`chemicalelements`](@ref) -> [`chemicalformula`](@ref).

# Keyword arguments
* `delim::Union{String, Char}` assigns the delimiter between each element.
* `unique::Bool` determines whether combines the elements to become unique or not when constructing formula from attribute [`chemicalelements`](@ref). 
It defaults to false for type [`Chemical`](@ref) and [`FormulaChemical`](@ref).
* `ischemical::Bool` determines whether the chemical is a chemical or a scheme. 
* `loss::Bool` determines whether the chemical is part of chemical loss, and signs are factored out from elements. 
"""
function chemicalformula(chemical::AbstractChemical; shallow = false, loss = false, ischemical = true, kwargs...) 
    result = getchemicalproperty(chemical, :formula, "")
    if isempty(result) && !shallow
        chemicalformula(chemicalelements(chemical; loss = false, shallow = true); loss, ischemical, kwargs...)
    else 
        reverse_formula(result, ischemical, loss)
    end
end

function chemicalformula(chemical::AbstractScheme; shallow = false, loss = false, ischemical = false, kwargs...) 
    result = getchemicalproperty(chemical, :formula, "")
    if isempty(result) && !shallow
        chemicalformula(chemicalelements(chemical; loss = false, shallow = true); loss, ischemical, kwargs...)
    else 
        reverse_formula(result, ischemical, loss)
    end
end

function chemicalformula(adduct_ion::AbstractAdductIon; loss = false, kwargs...)
    el = copy(chemicalelements(ioncore(adduct_ion); loss, kwargs...))
    nm = ncore(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    chemicalformula(gain_elements!(dictionary_elements(Dictionary, el), chemicalelements(ionadduct(adduct_ion); loss, kwargs...)); loss = !loss, kwargs...)
end

chemicalformula(chemical::AbstractChemicalWrapper; loss = false, kwargs...) = chemicalformula(chemical.chemical; loss, kwargs...)
chemicalformula(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`chemicalformula` cannot be defined for structural scheme."))

"""
    chemicalelements(chemical::AbstractChemicalsSchema; loss = false, kwargs...) -> Vector{Pair{String, Int}}

The elements of chemical entity of `chemical`.

# Generic Methods
* `AbstractChemical`: property search.
* `AbstractAdductIon`: elements combining core chemical and adduct.
* `AbstractChemicalWrapper`: chemical elements of field `chemical`.
* `AbstractScheme`: property search; change of elements.
* `AbstractStructuralScheme`: throw error.

# Specific Methods 
* Entity Level.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref) -> `reverse_elements`.
    * property: `:elements`.
    * default: `Pair{String, Int}[]`.
2. Function: [`chemicalformula`](@ref) -> [`chemicalelements`](@ref).

# Keyword Arguments 
* `loss::Bool` determines whether the chemical is part of chemical loss, and sign flips are propagated into elements.
"""
function chemicalelements(chemical::AbstractChemical; shallow = false, loss = false, kwargs...) 
    result = getchemicalproperty(chemical, :elements, Pair{String, Int}[])
    if isempty(result) && !shallow
        chemicalelements(chemicalformula(chemical; loss = false, shallow = true); loss, kwargs...)
    else 
        reverse_elements(result, loss)
    end
end

function chemicalelements(chemical::AbstractScheme; shallow = false, loss = false, kwargs...) 
    result = getchemicalproperty(chemical, :elements, Pair{String, Int}[])
    if isempty(result) && !shallow
        chemicalelements(chemicalformula(chemical; loss = false, ischemical = false, shallow = true); loss, kwargs...)
    else 
        reverse_elements(result, loss) 
    end
end

function chemicalelements(adduct_ion::AbstractAdductIon; loss = false, kwargs...) 
    el = copy(chemicalelements(ioncore(adduct_ion); loss, kwargs...))
    nm = ncore(adduct_ion)
    if nm > 1
        for id in eachindex(el)
            el[id] = first(el[id]) => nm * last(el[id])
        end
    end
    vcat(el, chemicalelements(ionadduct(adduct_ion); loss, kwargs...))
end

chemicalelements(chemical::AbstractChemicalWrapper; loss = false, kwargs...) = chemicalelements(chemical.chemical; loss, kwargs...)
chemicalelements(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`chemicalelements` cannot be defined for structural scheme."))

"""
    chemicalabbr(chemical::AbstractChemicalsSchema; verbose = true, n = 1, loss = false, bracket = true, kwargs...) -> String

The abbreviation of `chemical`. 

# Generic Methods
* `AbstractChemicalsSchema`: property search.
* `AbstractAdductIon`: abbreviation composed of abbreviation of core chemical and adduct.
* `AbstractChemicalWrapper`: abbreviation of field `chemical`.

# Specific Methods 
* Species/Transition Level.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref).
    * property: `:abbreviation`.
    * default: [`chemicalname`](@ref).

# Keyword arguments
* `verbose::Bool` determines whether includes all abbreviations or not for chemical species. If `verbose` is false, only the first (most abundant) chemical is included.
* `n::Int` determines number of chemical.
* `loss::Bool` determines whether the chemical is part of chemical loss, and sign flips are propagated into this chemical.
* `bracket::Bool` determines whether wrapping the output string with "[...]" and charge.
"""
function chemicalabbr(chemical::AbstractChemicalsSchema; n = 1, kwargs...) 
    abbr = getchemicalproperty(chemical, :abbreviation, "")
    isempty(abbr) && return chemicalname(chemical; n, kwargs...)
    string(n > 1 ? n : "", abbr)
end
chemicalabbr(chemical::AbstractChemicalWrapper; n = 1, kwargs...) = chemicalabbr(chemical.chemical; n, kwargs...)
function chemicalabbr(adduct_ion::AbstractAdductIon; loss = false, bracket = true, n = 1, kwargs...) 
    r = chemicalabbr(ioncore(adduct_ion); loss = bracket ? false : loss, bracket = false, n, kwargs...)
    a = chemicalabbr(ionadduct(adduct_ion); loss = bracket ? false : loss, bracket = false, n, kwargs...)
    n = ncore(adduct_ion)
    if occursin(" ", r)
        r = string("(", r, ")")
    end
    if n > 1 
        r = string(n, r)
    end
    bracket ? string("[", r, a, "]", charge_repr(charge(adduct_ion))) : string(r, a)
end

"""
    chemicalsmiles(chemical::AbstractChemical; kwargs...) -> String

The SMILES of chemical entity of `chemical`.

# Generic Methods
* `AbstractChemical`: property search.
* `AbstractAdductIon`: SMILES of the core chemical.
* `AbstractChemicalWrapper`: SMILES of field `chemical`.

# Specific Methods 
* Entity Level.
* `Isotopomers`: SMILES of the parent chemical.
* `Groupedisotopomers`: SMILES of the most abundant parent chemical.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref).
    * property: `:SMILES`.
    * default: `""`.
"""
chemicalsmiles(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :SMILES, "")
chemicalsmiles(adduct_ion::AbstractAdductIon; kwargs...) = chemicalsmiles(ioncore(adduct_ion); kwargs...) 
chemicalsmiles(chemical::AbstractChemicalWrapper; kwargs...) = chemicalsmiles(chemical.chemical; kwargs...)

"""
    ioncore(adduct_ion::AbstractAdductIon{S, T}; kwargs...) -> S

The core chemical of `adduct_ion`. 

# Generic Methods
* `AbstractAdductIon`: property search.

# Specific Methods 
* Entity Level.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref).
    * property: `:core`.
    * default: `TypeError`.
"""
ioncore(adduct_ion::AbstractAdductIon{S, T}; kwargs...) where {S, T} = getchemicalproperty(adduct_ion, :core, nothing)::S

"""
    ionadduct(adduct_ion::AbstractAdductIon{S, T}; kwargs...) -> T

The adduct of `adduct_ion`. 

# Generic Methods
* `AbstractAdductIon`: property search.

# Specific Methods 
* Entity Level.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref).
    * property: `:adduct`.
    * default: `TypeError`.
"""
ionadduct(adduct_ion::AbstractAdductIon{S, T}; kwargs...) where {S, T} = getchemicalproperty(adduct_ion, :adduct, nothing)::T

"""
    ncore(adduct_ion::AbstractAdductIon; kwargs...) -> Int

The number of core chemical. For instance, 2 for "[2M+H]+". 

# Generic Methods
* `AbstractAdductIon`: property search.

# Specific Methods 
* Entity Level.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref).
    * property: `:ncore`.
    * default: `1`.
"""
ncore(adduct_ion::AbstractAdductIon; kwargs...) = getchemicalproperty(adduct_ion, :ncore, 1)

"""
    charge(chemical::AbstractChemicalsSchema; loss = false, kwargs...) -> Int

The charge state of `chemical`; positive for cation and negative for anion. For instance, -2 for dianion, +3 for trication. 

# Generic Methods
* `AbstractChemicals`: property search.
* `AbstractAdductIon`: charge of core chemical and adduct.
* `AbstractChemicalWrapper`: charge of field `chemical`.
* `AbstractScheme`: property search; change of charges.
* `AbstractCompleteScheme`: change of charges of elemental scheme.
* `AbstractStructuralScheme`: throw error.

# Specific Methods 
* Entity Level.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref) -> sign flip if `loss`.
    * property: `:charge`.
    * default: `0`.

# Keword Arguments
* `loss::Bool` determines whether the chemical is part of chemical loss, and sign flips are factored out from elements.
"""
charge(chemical::AbstractChemicalsSchema; loss = false, kwargs...) = loss ? -getchemicalproperty(chemical, :charge, 0) : getchemicalproperty(chemical, :charge, 0)
charge(adduct_ion::AbstractAdductIon; loss = false, kwargs...) = ncore(adduct_ion) * charge(ioncore(adduct_ion); loss, kwargs...) + charge(ionadduct(adduct_ion); loss)
charge(chemical::AbstractChemicalWrapper; loss = false, kwargs...) = charge(chemical.chemical; loss, kwargs...)
charge(sch::AbstractCompleteScheme; loss = false, kwargs...) = charge(elementalscheme(sch); loss, kwargs...)
charge(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`charge` cannot be defined for structural scheme."))

"""
    ncharge(chemical::AbstractChemicalsSchema; kwargs...) -> Int

The number of charges of `chemical`.

# Generic Methods
* `AbstractChemicalsSchema`: absolute value of charge.

# Specific Methods 
* Entity Level.
"""
ncharge(chemical::AbstractChemicalsSchema; kwargs...) = abs(charge(chemical; kwargs...))

"""
    retentiontime(chemical::AbstractChemical; kwargs...) -> AbstractFloat

The retention time of `chemical`.

# Generic Methods
* `AbstractChemical`: property search.
* `AbstractChemicalWrapper`: retention time of field `chemical`.

# Specific Methods 
* Entity Level.
* `Isobars`: weighted mean of retention times of each chemical.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref).
    * property: `:retentiontime`.
    * default: `NaN`.
"""
retentiontime(chemical::AbstractChemical; kwargs...) = getchemicalproperty(chemical, :retentiontime, NaN)
retentiontime(chemical::AbstractChemicalWrapper; kwargs...) = retentiontime(chemical.chemical; kwargs...)

"""
    chemicalentity(chemical::AbstractChemical; kwargs...) -> AbstractChemical

The single chemical entity (having a single formula) from a chemical entity (i.e. itself), or chemical species. 

Attribute with Specific Methods marked as `Entity` indicates the `chemicalentity` is applied in the function for chemical type that is not originally a single entity. 

# Generic Methods
* `AbstractChemical`: itself.
* `AbstractChemicalWrapper`: wrapped chemical entity of field `chemical`.

# Specific Methods
* `Isobars`: `chemicalentity(first(chemical.chemicals))`, i.e. the most abundant entity.
* `ChemicalTransition`: the very beginning precursor.
* `Groupedisotopomers`: the most abundant isotopomer.
"""
chemicalentity(chemical::AbstractChemical; kwargs...) = chemical
chemicalentity(chemical::T; kwargs...) where {T<:AbstractChemicalWrapper} = T.name.wrapper(chemicalentity(chemical.chemical; kwargs...))

"""
    elementalscheme(scheme::AbstractScheme; kwargs...) -> Union{ElementalSchema, AbstractChemical}

The elemental scheme of `scheme`. 

Attribute with Specific Methods marked as `Entity` indicates the `elementalscheme` is applied in the function for scheme type that is not originally elemental. 

# Generic Methods
* `AbstractScheme`: property search.
* `AbstractStructuralScheme`: throw error.

# Specific Methods
* `ChemicalSchema`: elemental schema of all schema.
* `IsotopomerizedSchema`: elemental schema of all schema.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref).
    * property: `:elementalscheme`.
    * default: itself.
"""
elementalscheme(scheme::AbstractScheme; kwargs...) = getchemicalproperty(scheme, :elementalscheme, scheme, Union{ElementalSchema, AbstractChemical})
elementalscheme(scheme::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("Structural scheme requires specific mapping to elemental scheme."))

"""
    structuralscheme(scheme::AbstractScheme; kwargs...) -> AbstractScheme

The structural scheme of `scheme`. All elemental schema are regarded as structural schema as well.

# Generic Methods
* `AbstractScheme`: property search.

# Specific Methods
* `ChemicalSchema`: structural schema of all schema.
* `IsotopomerizedSchema`: structural schema of all schema.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref).
    * property: `:structuralscheme`.
    * default: itself.
"""
structuralscheme(scheme::AbstractScheme; kwargs...) = getchemicalproperty(scheme, :structuralscheme, scheme, AbstractScheme)

"""
    chemicalspecies(chemical::AbstractChemical; kwargs...) -> Vector{<:AbstractChemical}

The chemical species (represented by a vector) from a chemical entity or chemical species. 

Attribute with Specific Methods marked as `Species` indicates that the return value may represent attribute of species. 

# Generic Methods
* `AbstractChemical`: `[chemical]`.

# Specific Methods
* `Isobars`: `chemical.chemicals`.
* `ChemicalTransition{<:Isobars}`: a vector of chemical transition.
"""
chemicalspecies(chemical::AbstractChemical; kwargs...) = [chemical]

"""
    chemicaltransition(chemical::AbstractChemical; kwargs...) -> Vector{<:AbstractChemical}

The chemical entities that are analyzed in each stage of instrumental analysis.

Attribute with Specific Methods marked as `Transition` indicates the return value may represent attribute of transition. 

# Generic Methods
* `AbstractChemical`: `[chemical]`.

# Specific Methods
* `ChemicalTransition`: a vector of chemical entities.
* `Isobars`: the most abundant chemical transition.
"""
chemicaltransition(chemical::AbstractChemical; kwargs...) = [chemical]

"""
    chemicalparent(chemical::AbstractChemicalsSchema; kwargs...) -> AbstractChemical

The parent chemical without delocalized isotopes replacement.

# Generic Methods
* `AbstractChemical`: itself.
* `AbstractChemicalWrapper`: wrapped parent chemical of field `chemical`.
* `AbstractScheme`: itself.
* `AbstractStructuralScheme`: throw error.

# Specific Methods
* Entity Level.
* `ChemicalTransition`: [`ChemicalTransition`](@ref) of parent chemicals of each transition.
* `Groupedisotopmers`: field `parent`.
* `ElementalScheme`: scheme of chemical parent of chemical entity.
* `ChemicalSchema`: schema of chemical parent(s) of chemical entity(s).
* `IsotopomerizedSchema`: field `parent`.
* `Groupedisotopmerizedschema`: field `parent`.
"""
chemicalparent(cc::AbstractChemical; kwargs...) = cc
chemicalparent(chemical::T; kwargs...) where {T<:AbstractChemicalWrapper} = T.name.wrapper(chemicalparent(chemical.chemical; kwargs...))
chemicalparent(sch::AbstractScheme; kwargs...) = sch
chemicalparent(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`chemicalparent` cannot be defined for structural scheme."))

"""
    isotopomersisotopes(chemical::AbstractChemicalsSchema; loss = false, kwargs...) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of isotopomers.

# Generic Methods
* `AbstractChemical`: property search.
* `AbstractChemicalWrapper`: delocalized isotopes replacement of field `chemical`.
* `AbstractScheme`: `Pair{String, Int}[]`.
* `AbstractStructuralScheme`: throw error.
* `AbstractCompleteScheme`: delocalized isotopes replacement of elemental scheme.

# Specific Methods 
* Entity Level.
* `Isotopomers`: field `isotopes`.
* `Groupedisotopomers`: isotopes replacement of the most abundant isotopomers.
* `ElementalScheme`: change of delocalized isotopes replacement.
* `ChemicalSchema`: all changes of delocalized isotopes replacement.
* `IsotopomerizedSchema`: field `isotopes`.
* `Groupedisotopomerizedschema`: isotopes replacement of the most abundant isotopomers.

# Property Search Workflow
1. Function: [`getchemicalproperty`](@ref).
    * property: `:isotopomersisotopes`.
    * default: `Pair{String, Int}[]`.

# Keyword arguments
* `loss::Bool` determines whether the chemical is part of chemical loss, and sign flips are propagated into elements.
"""
isotopomersisotopes(cc::AbstractChemicalsSchema; loss = false, kwargs...) = reverse_elements(getchemicalproperty(cc, :isotopomersisotopes, Pair{String, Int}[]), loss)
isotopomersisotopes(cc::AbstractChemicalWrapper; loss = false, kwargs...) = isotopomersisotopes(cc.chemical; loss, kwargs...)
isotopomersisotopes(sch::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`isotopomersisotopes` cannot be defined for structural scheme."))
isotopomersisotopes(sch::AbstractCompleteScheme; loss = false, kwargs...) = isotopomersisotopes(elementalscheme(sch); loss, kwargs...) 

"""
    isotopomerstate(cc::AbstractChemical; isotope = "[13C]", ischemical = true, loss = false, kwargs...) -> Int 
    isotopomerstate(cc::AbstractScheme; isotope = "[13C]", ischemical = false, loss = false, kwargs...) -> Int 

The isotopomers state, i.e. equivalent number of `isotope`. 

# Generic Methods
* `AbstractChemicalsSchema`: isotopomers state calculated from `isotopomersisotopes`.

# Specific Methods
* Entity Level.

# Keyword arguments
* `ischemical::Bool` determines whether the chemical is a chemical or a scheme. 
* `loss::Bool` determines whether the chemical is part of chemical loss, and sign flips are propagated into `isotopomersisotopes`.
"""
isotopomerstate(cc::AbstractChemical; isotope_unit = nothing, isotope = "[13C]", loss = false, ischemical = true, kwargs...) = _isotopomerstate(isotopomersisotopes(cc), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss, ischemical, kwargs...)
isotopomerstate(cc::AbstractScheme; isotope_unit = nothing, isotope = "[13C]", loss = false, ischemical = false, kwargs...) = _isotopomerstate(isotopomersisotopes(cc), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss, ischemical, kwargs...)

"""
    groupedisotopomersisotopes(chemical::AbstractChemicalsSchema; loss = false, kwargs...) -> Vector{Vector{Pair{String, Int}}}

The delocalized isotopes replacement of each isotopomers.

# Generic Methods
* `AbstractChemical`: single-element vector of `isotopomersisotopes`.
* `AbstractScheme`: `groupedisotopomersisotopes` of elemental scheme.
* `AbstractStructuralScheme`: throw error.

# Specific Methods 
* Entity Level.
* `Groupedisotopomers`: isotopes replacements of each isotopomers.
* `ElementalScheme`: single-element vector of `isotopomersisotopes`.
* `ChemicalSchema`: `Pair{String, Int}[]`.
* `IsotopomerizedSchema`: single-element vector of `isotopomersisotopes`.
* `Groupedisotopomerizedschema`: isotopes replacements of each isotopomers.

# Keyword arguments
* `loss::Bool` determines whether the chemical is part of chemical loss, and sign flips are propagated into elements.
"""
groupedisotopomersisotopes(x::AbstractChemical; loss = false, kwargs...) = [isotopomersisotopes(x; loss, kwargs...)]
groupedisotopomersisotopes(x::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`groupedisotopomersisotopes` cannot be defined for structural scheme."))
groupedisotopomersisotopes(x::AbstractScheme; loss = false, kwargs...) = groupedisotopomersisotopes(elementalscheme(x); loss, kwargs...)

"""
    groupedisotopomersabundance(chemical::AbstractChemicalsSchema; loss = false, kwargs...) -> Vector{<:AbstractFloat}

The abundance of each isotopomers.

# Generic Methods
* `AbstractChemical`: `[1.0]`.
* `AbstractScheme`: `groupedisotopomersabundance` of elemental scheme.
* `AbstractStructuralScheme`: throw error.

# Specific Methods 
* Entity Level.
* `Groupedisotopomers`: field `abundance`.
* `ElementalScheme`: `groupedisotopomersabundance` of field `chemical`.
* `ChemicalSchema`: `[1.0]`.
* `IsotopomerizedSchema`: `[1.0]`.
* `Groupedisotopomerizedschema`: field `abundance`.
"""
groupedisotopomersabundance(x::AbstractChemical; kwargs...) = [1.0]
groupedisotopomersabundance(x::AbstractStructuralScheme; kwargs...) = throw(ArgumentError("`groupedisotopomersabundance` cannot be defined for structural scheme."))
groupedisotopomersabundance(x::AbstractScheme; kwargs...) = groupedisotopomersabundance(elementalscheme(x); kwargs...)
# truly formed chemical
"""
    analyzedchemical(chemical::AbstractChemicalsSchema; kwargs...) -> AbstractChemical

The single chemical entity that is directly analyzed at the very beginning of instrumental analysis.

# Generic Methods
* `AbstractChemical`: itself.
* `AbstractScheme`: throw error.

# Specific Methods
* Species Level.
* `ChemicalTransition`: the very begining precursor, but cannot be chemical gain or loss.
* `Isobars`: the most abundant analyzed chemical.
"""
analyzedchemical(cc::AbstractChemical; kwargs...) = cc
analyzedchemical(cc::AbstractScheme; kwargs...) = throw(ArgumentError("Scheme cannot be analyzed directly."))

"""
    detectedchemical(chemical::AbstractChemicalsSchema; kwargs...) -> AbstractChemical

The single chemical entity that is directly detected at the very end of instrumental analysis.

# Generic Methods
* `AbstractChemical`: itself or precursor with scheme (scheme is derived from `detectedchemical` of `precursor` and itself as product using [`completescheme`](@ref)).
* `AbstractScheme`: precursor with scheme. It requires a keyword argument `precursor`.

# Specific Methods
* Species Level.
* `ChemicalTransition`: the very ending product, but takes chemical gain and loss in consideration.
* `Isobars`: the most abundant detected chemical.
"""
function detectedchemical(cc::AbstractChemical; precursor = nothing, kwargs...)
    isnothing(precursor) && return cc 
    precursor = detectedchemical(precursor; kwargs...)
    detectedchemical(precursor, cc)
end
function detectedchemical(sch::AbstractScheme; precursor = nothing, kwargs...) 
    isnothing(precursor) && throw(ArgumentError("Scheme cannot be directly detected without precursor."))
    precursor = detectedchemical(precursor; kwargs...)
    detectedchemical(precursor, sch)
end
"""
    detectedisotopes(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{Pair{String, Int}}

The delocalized isotopes replacement of detected chemical. See [`detectedchemical`](@ref) for details.
"""
detectedisotopes(cc::AbstractChemicalsSchema; kwargs...) = isotopomersisotopes(detectedchemical(cc); kwargs...)
function detectedisotopes(sch::AbstractScheme; precursor = nothing, precursorisotopes = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorisotopes) && throw(ArgumentError("Scheme cannot be directly detected without precursor information."))
    pre = isnothing(precursorisotopes) ? detectedisotopes(precursor; kwargs...) : precursorisotopes
    gain_elements(pre, isotopomersisotopes(sch; kwargs...))
end

"""
    detectedcharge(chemical::AbstractChemicalsSchema; kwargs...) -> Int

The charge state of detected chemical. See [`detectedchemical`](@ref) for details.
"""
detectedcharge(cc::AbstractChemicalsSchema; kwargs...) = charge(detectedchemical(cc); kwargs...)
function detectedcharge(sch::AbstractScheme; precursor = nothing, precursorcharge = nothing, kwargs...) 
    isnothing(precursor) && isnothing(precursorcharge) && throw(ArgumentError("Scheme cannot be directly detected without precursor information."))
    (isnothing(precursorcharge) ? detectedcharge(precursor; kwargs...) : precursorcharge) + charge(sch; kwargs...)
end

# MS representation of chemical
"""
    inputchemical(chemical::AbstractChemicalsSchema; kwargs...) -> AbstractChemical

The single chemical entity that is the input at the very beginning of instrumental analysis. 

It is equivalent to [`analyzedchemical`](@ref), except that schema are accepted. 
"""
inputchemical(cc::AbstractChemicalsSchema; kwargs...) = cc

"""
    outputchemical(chemical::AbstractChemicalsSchema; kwargs...) -> AbstractChemicalsSchema

The single chemical entity that is the output of the very ending of instrumental analysis.

It is equivalent to [`detectedchemical`](@ref), except that schema are kept unchanged. 
"""
outputchemical(cc::AbstractChemicalsSchema; kwargs...) = cc

"""
    seriesanalyzedchemical(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{<:AbstractChemical}

The chemical entities that are directly analyzed in each stage of instrumental analysis. 

It is equivalent to [`chemicaltransition`](@ref), except that schema are transformed to detected chemicals.

# Generic Methods
* `AbstractChemicalsSchema`: `[chemical]`.

# Specific Methods
* Transition Level.
* `ChemicalTransition`: a vector of analyzed chemicals. See [`analyzedchemical`](@ref) for details.
"""
seriesanalyzedchemical(cc::AbstractChemicalsSchema; kwargs...) = [cc]

"""
    seriesanalyzedisotopes(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{Vector{Pair{String, Int}}}

The delocalized isotopes replacement of serially analyzed chemical. See [`seriesanalyzedchemical`](@ref) for details.
"""
seriesanalyzedisotopes(cc::AbstractChemicalsSchema; kwargs...) = [isotopomersisotopes(c; kwargs...) for c in seriesanalyzedchemical(cc)]

"""
    seriesanalyzedcharge(chemical::AbstractChemicalsSchema; kwargs...) -> Vector{Int}

The charge states of serially analyzed chemical. See [`seriesanalyzedchemical`](@ref) for details.
"""
seriesanalyzedcharge(cc::AbstractChemicalsSchema; kwargs...) = [charge(c; kwargs...) for c in seriesanalyzedchemical(cc)]

"""
    msstage(chemical::AbstractChemical; kwargs...) -> Int

Number of stages of MS the chemical has been through.

# Generic Methods
* `AbstractChemical`: `1`.

# Specific Methods
* `ChemicalTransition`: `length(chemical.transition)`.
"""
msstage(cc::AbstractChemical; kwargs...) = 1

"""
    mmi(chemical::AbstractChemicalsScheme; loss = false, kwargs...) -> AbstractFloat

The monoisotopic mass of `chemical`.

# Generic Methods
* `AbstractChemical`: calculated from [`chemicalelements`](@ref) and [`charge`](@ref).
* `AbstractScheme`: calculated from [`chemicalelements`](@ref) and [`charge`](@ref); change of monoisotopic mass.
* `AbstractCompleteScheme`: change of monoisotopic mass of elemental scheme.

# Specific Methods
* Entity Level.
* `Isobars`: weighted mean of monoisotopic masses.
* `Groupedisotopomers`: weighted mean of monoisotopic masses.
* `ChemicalTransition`: monoisotopic mass of `analyzedchemical`.
* `Groupedisotopomerizedschema`: weighted mean of change of monoisotopic masses.

# Keword Arguments
* `loss::Bool` determines whether the chemical is part of chemical loss, and sign flips are factored out from elements.
"""
mmi(cc::AbstractChemicalsSchema; loss = false) = mmi(chemicalelements(cc), charge(cc); loss)
mmi(sch::AbstractCompleteScheme; loss = false, kwargs...) = mmi(elementalscheme(sch); loss, kwargs...)

"""
    molarmass(chemical::AbstractChemicalsSchema; loss = false, kwargs...) -> AbstractFloat

The molar mass of `chemical`.

# Generic Methods
* `AbstractChemical`: calculated from [`chemicalelements`](@ref) and [`charge`](@ref).
* `AbstractScheme`: calculated from [`chemicalelements`](@ref) and [`charge`](@ref); change of molar mass.
* `AbstractCompleteScheme`: change of molar mass of elemental scheme.

# Specific Methods
* Entity Level.
* `Isobars`: weighted mean of molar masses.
* `Groupedisotopomers`: weighted mean of molar masses.
* `ChemicalTransition`: molar mass of `analyzedchemical`.
* `Groupedisotopomerizedschema`: weighted mean of change of molar masses.

# Keword Arguments
* `loss::Bool` determines whether the chemical is part of chemical loss, and sign flips are factored out from elements.
"""
molarmass(cc::AbstractChemicalsSchema; loss = false) = molarmass(chemicalelements(cc), charge(cc); loss)
molarmass(sch::AbstractCompleteScheme; loss = false, kwargs...) = molarmass(elementalscheme(sch); loss, kwargs...)

"""
    mz(chemical::AbstractChemical[, adduct]; loss = false, kwargs...) -> AbstractFloat
    mz(chemical::AbstractScheme; loss = false, kwargs...) -> AbstractFloat

The mass to charge ratio (m/z) of charged chemical or chemical with adduct. It is equivalent to `mmi(charged_chemical) / ncharge(charged_chemical)`.

# Generic Methods
* `AbstractChemical`: calculated from [`mmi`](@ref) and [`charge`](@ref).
* `AbstractAdductIon`: calculated from [`mmi`](@ref) and [`charge`](@ref) of both core chemical and adduct.
* `AbstractScheme`: calculated from [`mmi`](@ref) and [`charge`](@ref); change of m/z.
* `AbstractCompleteScheme`: change of m/z of elemental scheme.

# Specific Methods
* Entity Level.
* `Isobars`: weighted mean of m/z.
* `Groupedisotopomers`: weighted mean of m/z.
* `ChemicalTransition`: m/z of `analyzedchemical`.
* `Groupedisotopomerizedschema`: weighted mean of change of m/z.

# Keword Arguments
* `loss:Bool` determines whether the chemical is part of chemical loss, and sign flips are factored out from elements.
"""
mz(charged_cc::AbstractChemicalsSchema; loss = false, kwargs...) = charge(charged_cc) == 0 ? NaN : mmi(charged_cc; loss, kwargs...) / ncharge(charged_cc)
mz(cc::AbstractChemical, adduct; loss = false, kwargs...) = mz(ionize(cc; parse_adduct(adduct)...); loss, kwargs...)
function mz(adduct_ion::AbstractAdductIon; loss = false, kwargs...) 
    (ncore(adduct_ion) * mmi(ioncore(adduct_ion); loss, kwargs...) + mmi(ionadduct(adduct_ion); loss, kwargs...)) / ncharge(adduct_ion)
end
mz(sch::AbstractCompleteScheme; loss = false, kwargs...) = mz(elementalscheme(sch); loss, kwargs...)


# list all propertys
# add new propertys
# function chemicalpropertynames() end