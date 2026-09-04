include(joinpath("elements", "consts.jl"))
include(joinpath("elements", "property.jl"))
include(joinpath("elements", "manipulation.jl"))

"""
    set_element!(element, mass, abundance; minor_name = nothing)

Update or insert `element`. 

# Arguments
* `element::AbstractString`: element name.
* `mass::Vector`: atomic mass of all isotopes.
* `abundance::Vector`: natural abundance of all isotopes.
* `minor_name`: customized minor element names.
"""
function set_element!(element::AbstractString, mass, abundance; minor_name = nothing)
    isotopes = isnothing(minor_name) ? map(mass) do m 
        n = round(Int, m)
        string("[", n, element, "]")
    end : minor_name
    _, i = findmax(abundance)
    elements_parents()[element] = element
    elements_decodes()[element] = element
    elements_mass()[element] = mass[i]
    elements_abundance()[element] = abundance[i]
    ee = encode_isotopes.(isotopes)
    for (e, i, m, a) in zip(ee, isotopes, mass, abundance) 
        elements_decodes()[e] = i 
        elements_mass()[i] = m 
        elements_abundance()[i] = a 
        elements_parents()[i] = element
    end
    id = sortperm(abundance; rev = true)
    elements_isotopes()[element] = isotopes[id]
    MASS
end

"""
    chemicalformula(elements::Dict; delim = "", unique = true, ischemical = true, loss = false) -> String
    chemicalformula(elements::Dictionary; delim = "", unique = true, ischemical = true, loss = false) -> String
    chemicalformula(elements::Vector{<:Pair}; delim = "", unique = true, ischemical = true, loss = false) -> String
    chemicalformula(elements::ElementsVector; delim = "", unique = true, ischemical = true, loss = false) -> String

Create chemical formula using given element-number pairs. 

# Arguments
* `delim::Union{String, Char}` assigns the delimiter between each element.
* `unique::Bool` determines whether combines the elements to become unique or not.
* `ischemical::Bool` determines whether the chemical is a chemical or a scheme. 
* `loss::Bool` determines whether the chemical is part of chemical loss, and signs are factored out from elements. 
"""
function chemicalformula(elements::Vector{<:Pair}; delim = "", unique = true, loss = false, ischemical = true)
    if unique 
        chemicalformula(dictionary_elements(Dictionary, elements); unique = false, loss, ischemical, delim)
    elseif all(x -> last(x) >= 0, elements)
        string(ischemical ? "" : loss ? "-" : "+", join((v == 1 ? k : string(k, v) for (k, v) in elements if v != 0), delim))
    elseif all(x -> last(x) <= 0, elements)
        string(loss ? "+" : "-", join((v == -1 ? k : string(k, abs(v)) for (k, v) in elements if v != 0), delim))
    elseif ischemical
        chemicalformula(dictionary_elements(Dictionary, elements); unique = false, loss, ischemical, delim)
    else
        elements_pos = filter(x -> last(x) > 0, elements)
        elements_neg = filter(x -> last(x) < 0, elements)
        string(loss ? "+" : "-", join((v == -1 ? k : string(k, abs(v)) for (k, v) in elements_neg), delim), loss ? "-" : "+", join((v == 1 ? k : string(k, v) for (k, v) in elements_pos), delim))
    end
end

chemicalformula(elements::Dict; delim = "", unique = true, loss = false, ischemical = true) = _chemicalformula(elements; delim, unique, loss, ischemical)
chemicalformula(elements::Dictionary; delim = "", unique = true, loss = false, ischemical = true) = _chemicalformula(pairs(elements); delim, unique, loss, ischemical)
chemicalformula(elements::ElementsVector; delim = "", unique = true, loss = false, ischemical = true) = _chemicalformula(elements; delim, unique, loss, ischemical)

function _chemicalformula(elements; delim = "", unique = false, loss = false, ischemical = true)
    if all(x -> last(x) >= 0, elements)
        string(ischemical ? "" : loss ? "-" : "+", join((v == 1 ? k : string(k, v) for (k, v) in elements if v != 0), delim))
    elseif all(x -> last(x) <= 0, elements)
        string(loss ? "+" : "-", join((v == -1 ? k : string(k, abs(v)) for (k, v) in elements if v != 0), delim))
    else
        elements_pos = filter(x -> last(x) > 0, elements)
        elements_neg = filter(x -> last(x) < 0, elements)
        string(loss ? "+" : "-", join((v == -1 ? k : string(k, abs(v)) for (k, v) in elements_neg), delim), loss ? "-" : "+", join((v == 1 ? k : string(k, v) for (k, v) in elements_pos), delim))
    end
end

"""
    chemicalelements(formula::AbstractString; loss = false) -> Vector{Pair{String, Int}}

Create element-number pairs from chemical formula. 

# Arguments
* `loss::Bool` determines whether the chemical is part of chemical loss, and sign flips are propagated into elements.
"""
function chemicalelements(formula::AbstractString; loss = false, kwargs...)
    fs = split(formula, "+")
    v = Vector{Pair{String, Int}}[]
    for f in fs 
        fns = split(f, "-")
        fp = popfirst!(fns)
        isempty(fp) || push!(v, [elements_decodes()[k] => v * (loss ? -1 : 1) for (k, v) in parse_compound(encode_isotopes(fp))])
        for fn in fns
            isempty(fn) || push!(v, [elements_decodes()[k] => v * (loss ? 1 : -1) for (k, v) in parse_compound(encode_isotopes(fn))])
        end
    end
    vcat(v...)::Vector{Pair{String, Int}}
end

function encode_isotopes(formula::AbstractString)
    f = string(formula)
    f2 = f
    for i in eachmatch(r"\[(\d*)([^\]]*)\]", f)
        m, e = i
        delta = isempty(m) ? 0 : (parse(Int, m) - round(Int, elements_mass()[e]))
        e = delta > 0 ? string(e, "it") * "n" ^ delta :
            delta < 0 ? string(e, "it") * "p" ^ abs(delta) : string(e, "itz")
        f2 = replace(f2, i.match => e)
    end
    f2
end

chemicalformula(cc::Chemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalformula(cc::FormulaChemical; unique = false, kwargs...) = chemicalformula(cc.elements; unique, kwargs...)
chemicalformula(isobars::Isobars; kwargs...) = chemicalformula(chemicalentity(isobars); kwargs...)::String
function chemicalformula(x::Isotopomers; kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); loss = false))
    chemicalformula(isotopeelements(elements, x.isotopes); kwargs...)
end

function isotopeelements(elements, isotopes)
    for (k, v) in isotopes
        e = get(elements_parents(), k, k) 
        k == e && continue 
        v == 0 && continue
        elements[e] -= v 
        get!(elements, k, 0)
        elements[k] += v 
    end
    elements
end
isotopeelements_vec(elements, isotopes) = collect(pairs(isotopeelements(elements, isotopes)))

function chemicalformula(x::Groupedisotopomers; kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); loss = false))
    chemicalformula(isotopeelements(elements, x.isotopes[begin]); kwargs...)
end
chemicalformula(ct::ChemicalTransition; kwargs...) = chemicalformula(chemicalentity(ct); kwargs...)::String

chemicalformula(sch::AbstractCompleteScheme; kwargs...) = chemicalformula(elementalscheme(sch); kwargs...)
chemicalformula(sch::ElementalScheme{false}; loss = false, kwargs...) = chemicalformula(sch.chemical; loss = !loss, kwargs..., ischemical = false)
chemicalformula(sch::ElementalScheme{true}; loss = false, kwargs...) = chemicalformula(sch.chemical; loss, kwargs..., ischemical = false)
chemicalformula(x::ChemicalSchema; kwargs...) = chemicalformula(chemicalelements(x; loss = false); kwargs..., ischemical = false)
function chemicalformula(x::IsotopomerizedSchema; kwargs...)
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); loss = false))
    chemicalformula(isotopeelements(elements, x.isotopes); kwargs..., ischemical = false)
end
function chemicalformula(x::Groupedisotopomerizedschema; kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); loss = false))
    chemicalformula(isotopeelements(elements, x.isotopes[begin]); kwargs..., ischemical = false)
end

function reverse_formula(x, ischemical, loss) 
    if isempty(x)
        x
    elseif ischemical 
        x 
    elseif startswith(x, r"[^+-]")
        loss ? string("-", replace(x, "+" => "-", "-" => "+")) : string("+", x)
    else
        loss ? replace(x, "+" => "-", "-" => "+") : x
    end
end
    
reverse_elements(x::ElementsVector, loss) = loss ? ElementsVector(x.elements, [-v for v in x.numbers]) : x
reverse_elements(x::Vector{<:Pair}, loss) = loss ? [k => -v for (k, v) in x] : x
reverse_elements(x::Dict, loss) = loss ? Dict(k => -v for (k, v) in x) : x
reverse_elements(x::Dictionary, loss) = loss ? Dictionary(keys(x), [-v for v in x]) : x

chemicalelements(cc::Chemical; loss = false, kwargs...) = reverse_elements(cc.elements, loss)
chemicalelements(cc::FormulaChemical; loss = false, kwargs...) = reverse_elements(cc.elements, loss)
chemicalelements(isobars::Isobars; kwargs...) = chemicalelements(chemicalentity(isobars); kwargs...)::Vector{Pair{String, Int}}
function chemicalelements(x::Isotopomers; loss = false, kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs..., loss = false))
    reverse_elements(isotopeelements_vec(elements, x.isotopes), loss)
end

function chemicalelements(x::Groupedisotopomers; loss = false, kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs..., loss = false))
    reverse_elements(isotopeelements_vec(elements, x.isotopes[begin]), loss)
end
chemicalelements(ct::ChemicalTransition; loss = false, kwargs...) = chemicalelements(chemicalentity(ct); loss, kwargs...)::Vector{Pair{String, Int}}

chemicalelements(sch::AbstractCompleteScheme; kwargs...) = chemicalelements(elementalscheme(sch); kwargs...)
chemicalelements(sch::ElementalScheme{false}; loss = false, kwargs...) = chemicalelements(sch.chemical; loss = !loss, kwargs...) 
chemicalelements(sch::ElementalScheme{true}; loss = false, kwargs...) = chemicalelements(sch.chemical; loss, kwargs...) 
function chemicalelements(x::IsotopomerizedSchema; loss = false, kwargs...)
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs..., loss = false))
    reverse_elements(isotopeelements_vec(elements, x.isotopes), loss)
end
chemicalelements(x::ChemicalSchema; kwargs...) = vcat((repeat(chemicalelements(k; kwargs...), v) for (k, v) in zip(x.schema, x.number))...)
function chemicalelements(x::Groupedisotopomerizedschema; loss = false, kwargs...) 
    elements = dictionary_elements(Dictionary, chemicalelements(chemicalparent(x); kwargs..., loss = false))
    reverse_elements(isotopeelements_vec(elements, x.isotopes[begin]), loss)
end

isotopomersisotopes(isobars::Isobars; kwargs...) = isotopomersisotopes(chemicalentity(isobars); kwargs...)::Vector{Pair{String, Int}}
isotopomersisotopes(isotopomers::Isotopomers; loss = false, kwargs...) = collect(reverse_elements(isotopomers.isotopes, loss))
isotopomersisotopes(isotopomers::Groupedisotopomers; loss = false, kwargs...) = collect(reverse_elements(isotopomers.isotopes[begin], loss))
isotopomersisotopes(ct::ChemicalTransition; kwargs...) = isotopomersisotopes(chemicalentity(ct); kwargs...)::Vector{Pair{String, Int}}

isotopomersisotopes(sch::ElementalScheme{true}; loss = false, kwargs...) = isotopomersisotopes(sch.chemical; loss, kwargs...)
isotopomersisotopes(sch::ElementalScheme{false}; loss = false, kwargs...) = isotopomersisotopes(sch.chemical; loss = !loss, kwargs...)
isotopomersisotopes(x::IsotopomerizedSchema; loss = false, kwargs...) = collect(reverse_elements(x.isotopes, loss))
isotopomersisotopes(x::Groupedisotopomerizedschema; loss = false, kwargs...) = collect(reverse_elements(x.isotopes[begin], loss))

isotopomerstate(sch::ElementalScheme{true}; isotope_unit = nothing, isotope = "[13C]", loss = false, kwargs...) = _isotopomerstate(isotopomersisotopes(sch; loss = false), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss, kwargs..., ischemical = false)
isotopomerstate(sch::ElementalScheme{false}; isotope_unit = nothing, isotope = "[13C]", loss = false, kwargs...) = _isotopomerstate(isotopomersisotopes(sch; loss = false), isnothing(isotope_unit) ? elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]] : isotope_unit; loss = !loss, kwargs..., ischemical = false)

function _isotopomerstate(isotopes::Vector, isotope_unit; ischemical = true, loss = false)
    ds = 0
    if ischemical || !loss
        for (e, n) in isotopes
            ds += (elements_mass()[e] - elements_mass()[elements_parents()[e]]) * n
        end
    else
        for (e, n) in isotopes
            ds += (elements_mass()[e] - elements_mass()[elements_parents()[e]]) * (-n)
        end
    end
    round(Int, ds / isotope_unit)
end

groupedisotopomersisotopes(x::ElementalScheme{true}; loss = false, kwargs...) = groupedisotopomersisotopes(x.chemical; loss, kwargs...)
groupedisotopomersisotopes(x::ElementalScheme{false}; loss = false, kwargs...) = groupedisotopomersisotopes(x.chemical; loss = !loss, kwargs...)
groupedisotopomersisotopes(x::ChemicalSchema; loss = false, kwargs...) = Pair{String, Int}[]
groupedisotopomersisotopes(x::Groupedisotopomers; loss = false, kwargs...) = [collect(reverse_elements(y, loss)) for y in x.isotopes]
groupedisotopomersisotopes(x::Groupedisotopomerizedschema; loss = false, kwargs...) = [collect(reverse_elements(y, loss)) for y in x.isotopes]

groupedisotopomersabundance(x::ElementalScheme; kwargs...) = groupedisotopomersabundance(x.chemical; kwargs...)
groupedisotopomersabundance(x::ChemicalSchema; kwargs...) = [1.0]
groupedisotopomersabundance(x::Groupedisotopomers; kwargs...) = x.abundance
groupedisotopomersabundance(x::Groupedisotopomerizedschema; kwargs...) = x.abundance