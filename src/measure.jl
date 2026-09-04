"""
    mmi(formula::AbstractString, net_charge = 0; loss = false) -> AbstractFloat
    mmi(elements, net_charge = 0; loss = false) -> AbstractFloat

Monoisotopic mass of formula, and elements.
"""
function mmi(elements::Union{<:Vector{<:Pair}, <:Dict}, net_charge = 0; loss = false)
    # Vector of el => #el
    weight = elements_mass()[""]
    for (el, n) in elements
        weight += elements_mass()[el] * n
    end
    weight -= net_charge * ME
    loss ? -weight : weight
end

function mmi(formula::AbstractString, net_charge = 0; loss = false)
    # if any prefix number
    n = match(r"^[0-9]+", formula)
    if isnothing(n)
        mmi(chemicalelements(formula), net_charge; loss)
    else
        formula = replace(formula, n.match => "")
        mmi(chemicalelements(formula), net_charge; loss) * parse(Int, n.match)
    end
end

"""
    molarmass(formula::AbstractString, net_charge = 0; loss = false) -> AbstractFloat
    molarmass(elements, net_charge = 0; loss = false) -> AbstractFloat

Molar mass of formula, and elements.
"""
function molarmass(elements::Union{<:Vector{<:Pair}, <:Dict}, net_charge = 0; loss = false)
    # Vector of el => #el
    weight = elements_mass()[""]
    for (el, n) in elements
        if iselement(el)
            for i in elements_isotopes()[el]
                weight += elements_mass()[i] * n * elements_abundance()[i]
            end
        else
            weight += elements_mass()[el] * n
        end
    end
    weight -= net_charge * ME
    loss ? -weight : weight
end

function molarmass(formula::AbstractString, net_charge = 0; loss = false)
    # if any prefix number
    n = match(r"^[0-9]+", formula)
    if isnothing(n)
        molarmass(chemicalelements(formula), net_charge; loss)
    else
        formula = replace(formula, n.match => "")
        molarmass(chemicalelements(formula), net_charge; loss) * parse(Int, n.match)
    end
end

"""
    value_error(x, y)

Error function of difference; `y - x`.
"""
value_error(x, y) = y - x

"""
    relative_error(x, y)

Error function of relative difference (relative to true value x); `(y - x) / x`.
"""
relative_error(x, y) = (y - x) / x

"""
    relative_error_mean(x, y)

Error function of relative difference (relative to mean); `(y - x) / ((x + y) / 2)`.
"""
relative_error_mean(x, y) = 2 * (y - x) / (x + y) 

"""
    percentage_error(x, y)

Error function of relative difference in percentage (relative to true value x); `(y - x) / x * 100`.
"""
percentage_error(x, y) = (y - x) / x * 100

"""
    percentage_error_mean(x, y)

Error function of relative difference in percentage (relative to mean); `(y - x) / ((x + y) / 2) * 100`.
"""
percentage_error_mean(x, y) = 2 * (y - x) / (x + y) * 100

"""
    ppm_error(x, y)

Error function of relative difference in ppm (relative to true value x); `(y - x) / x * 1e6`.
"""
ppm_error(x, y) = (y - x) / x * 1e6

"""
    ppm_error_mean(x, y)

Error function of relative difference in ppm (relative to mean); `(y - x) / ((x + y) / 2) * 1e6`.
"""
ppm_error_mean(x, y) = 2 * (y - x) / (x + y) * 1e6

"""
    measure_name(fn[, error]) -> String 

Common name of measurement `fn` with or without `error`.
"""
measure_name(fn) = repr(fn)
measure_name(fn::typeof(retentiontime)) = "RT"
measure_name(fn::typeof(mz)) = "MZ"
measure_name(fn::typeof(mmi)) = "Mmi"
measure_name(fn::typeof(molarmass)) = "M"
measure_name(fn::typeof(charge)) = "Z"
measure_name(fn::typeof(ncharge)) = "|Z|"
measure_name(s::AbstractString) = string(s)
measure_name(s::Symbol) = string(s)
measure_name(fn, error::typeof(value_error)) = string("Δ", measure_name(fn))
measure_name(fn, error::typeof(relative_error)) = string("Δ", measure_name(fn), "/", measure_name(fn))
measure_name(fn, error::typeof(relative_error_mean)) = string("Δ", measure_name(fn), "/", measure_name(fn))
measure_name(fn, error::typeof(percentage_error)) = string("Δ", measure_name(fn), "/", measure_name(fn), "(%)")
measure_name(fn, error::typeof(percentage_error_mean)) = string("Δ", measure_name(fn), "/", measure_name(fn), "(%)")
measure_name(fn, error::typeof(ppm_error)) = string("Δ", measure_name(fn), "/", measure_name(fn), "(ppm)")
measure_name(fn, error::typeof(ppm_error_mean)) = string("Δ", measure_name(fn), "/", measure_name(fn), "(ppm)")

"""
    measure_error(fn) -> Vector{<:Function}

Default error functions for measurement `fn`.
"""
measure_error(::typeof(retentiontime)) = [value_error]
measure_error(fn) = [value_error]

charge(isobars::Isobars; kwargs...) = charge(chemicalentity(isobars); kwargs...)
charge(isotopomers::Isotopomers; kwargs...) = charge(chemicalparent(isotopomers); kwargs...)
charge(isotopomers::Groupedisotopomers; kwargs...) = charge(chemicalparent(isotopomers); kwargs...)
charge(ct::ChemicalTransition; kwargs...) = charge(chemicalentity(ct); kwargs...)::Int

charge(x::ElementalScheme{false}; loss = false, kwargs...) = charge(x.chemical; loss = !loss, kwargs...) 
charge(x::ElementalScheme{true}; loss = false, kwargs...) = charge(x.chemical; loss, kwargs...) 
charge(x::IsotopomerizedSchema; kwargs...) = charge(chemicalparent(x); kwargs...)
charge(x::ChemicalSchema; kwargs...) = sum(charge(k; kwargs...) * v for (k, v) in zip(x.schema, x.number))
charge(x::Groupedisotopomerizedschema; kwargs...) = charge(chemicalparent(x); kwargs...)

retentiontime(isobars::Isobars; kwargs...) = _isobar_species_attr(retentiontime, isobars; kwargs...)
retentiontime(isotopomers::Isotopomers; kwargs...) = retentiontime(chemicalparent(isotopomers); kwargs...)
retentiontime(isotopomers::Groupedisotopomers; kwargs...) = retentiontime(chemicalparent(isotopomers); kwargs...)
retentiontime(ct::ChemicalTransition; kwargs...) = retentiontime(chemicalentity(ct); kwargs...)

msstage(isobars::Isobars{<:ChemicalTransition}; kwargs...) = only(unique(msstage.(chemicalspecies(isobars); kwargs...)))
msstage(ct::ChemicalTransition; kwargs...) = length(ct.transition)

_isobar_species_attr(fn, isobars::Isobars, args...; kwargs...) = mean(fn.(chemicalspecies(isobars), args...; kwargs...), weights(isobars.abundance))

element_mmi(x) = elements_mass()[x]

function _mass_isotope(isotopes; loss = false, kwargs...) 
    m = sum((element_mmi(x) - element_mmi(elements_parents()[x])) * n for (x, n) in isotopes; init = 0)
    loss ? -m : m
end
function _mass_isotope(isotopes, abundance; loss = false, kwargs...) 
    m = mean([sum((element_mmi(x) - element_mmi(elements_parents()[x])) * n for (x, n) in iso; init = 0) for iso in isotopes], weights(abundance))
    loss ? -m : m
end

mmi(isobars::Isobars; kwargs...) = _isobar_species_attr(mmi, isobars; kwargs...)
mmi(x::Isotopomers; kwargs...) = mmi(chemicalparent(x); kwargs...) + _mass_isotope(x.isotopes; kwargs...)
mmi(x::Groupedisotopomers; kwargs...) = mmi(chemicalparent(x); kwargs...) + _mass_isotope(x.isotopes, x.abundance; kwargs...)
mmi(ct::ChemicalTransition; kwargs...) = mmi(analyzedchemical(ct); kwargs...)::float(Int)

mmi(x::ElementalScheme{false}; loss = false, kwargs...) = mmi(x.chemical; loss = !loss, kwargs...) 
mmi(x::ElementalScheme{true}; loss = false, kwargs...) = mmi(x.chemical; loss, kwargs...)
mmi(x::ChemicalSchema; kwargs...) = sum(mmi(k; kwargs...) * v for (k, v) in zip(x.schema, x.number)) 
mmi(x::IsotopomerizedSchema; kwargs...) = mmi(chemicalparent(x); kwargs...) + _mass_isotope(x.isotopes; kwargs...)
mmi(x::Groupedisotopomerizedschema; kwargs...) = mmi(chemicalparent(x); kwargs...) + _mass_isotope(x.isotopes, x.abundance; kwargs...)

vec_mmi_fix(x, y; kwargs...) = [mmi(m) + y for m in x]

molarmass(isobars::Isobars; kwargs...) = _isobar_species_attr(molarmass, isobars; kwargs...)
molarmass(x::Isotopomers; kwargs...) = molarmass(chemicalparent(x); kwargs...) + _mass_isotope(x.isotopes; kwargs...)
molarmass(x::Groupedisotopomers; kwargs...) = molarmass(chemicalparent(x); kwargs...) + _mass_isotope(x.isotopes, x.abundance; kwargs...)
molarmass(ct::ChemicalTransition; kwargs...) = molarmass(analyzedchemical(ct); kwargs...)::float(Int)

molarmass(x::ElementalScheme{false}; loss = false, kwargs...) = molarmass(x.chemical; loss = !loss, kwargs...) 
molarmass(x::ElementalScheme{true}; loss = false, kwargs...) = molarmass(x.chemical; loss, kwargs...) 
molarmass(x::ChemicalSchema; kwargs...) = sum(molarmass(k; kwargs...) * v for (k, v) in zip(x.schema, x.number))
molarmass(x::IsotopomerizedSchema; kwargs...) = molarmass(chemicalparent(x); kwargs...) + _mass_isotope(x.isotopes; kwargs...)
molarmass(x::Groupedisotopomerizedschema; kwargs...) = molarmass(chemicalparent(x); kwargs...) + _mass_isotope(x.isotopes, x.abundance; kwargs...)

mz(ct::ChemicalTransition; kwargs...) = mz(analyzedchemical(ct); kwargs...)::float(Int)
mz(ct::ChemicalTransition, adduct; kwargs...) = mz(analyzedchemical(ct), adduct; kwargs...)::float(Int)