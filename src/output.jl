"""
    defaultname(chemical::AbstractChemicalsSchema)

Default name of `chemical` if no attribute and specific method is not implemented for `chemicalname`.
"""
defaultname(::T) where {T<:AbstractChemical} = string("Chemical::", T)
defaultname(::T) where {T<:AbstractScheme} = string("Scheme::", T)

post_decorator(::ElementalScheme{false}; delim = "") = delim
post_decorator(::ElementalScheme{false, <:FormulaChemical}; delim = "") = ""
post_decorator(::ElementalScheme{false, <:Isotopomers{<:FormulaChemical}}; delim = "") = ""
post_decorator(::ElementalScheme{false, <:Isotopomers{<:AbstractAdductIon{<:FormulaChemical}}}; delim = "") = ""
post_decorator(::ElementalScheme{false, <:Groupedisotopomers{<:FormulaChemical}}; delim = "") = ""
post_decorator(::ElementalScheme{false, <:Groupedisotopomers{<:AbstractAdductIon{<:FormulaChemical}}}; delim = "") = ""
post_decorator(::ElementalScheme{true}; delim = "") = delim
post_decorator(::ElementalScheme{true, <:FormulaChemical}; delim = "") = ""
post_decorator(::ElementalScheme{true, <:Isotopomers{<:FormulaChemical}}; delim = "") = ""
post_decorator(::ElementalScheme{true, <:Isotopomers{<:AbstractAdductIon{<:FormulaChemical}}}; delim = "") = ""
post_decorator(::ElementalScheme{true, <:Groupedisotopomers{<:FormulaChemical}}; delim = "") = ""
post_decorator(::ElementalScheme{true, <:Groupedisotopomers{<:AbstractAdductIon{<:FormulaChemical}}}; delim = "") = ""

decorator(::ElementalScheme{false}; loss = false, delim = "") = string(delim, loss ? "Gain_" : "Loss_")
decorator(::ElementalScheme{false, <:FormulaChemical}; loss = false, delim = "") = loss ? "+" : "-"
decorator(::ElementalScheme{false, <:Isotopomers{<:FormulaChemical}}; loss = false, delim = "") = loss ? "+" : "-"
decorator(::ElementalScheme{false, <:Isotopomers{<:AbstractAdductIon{<:FormulaChemical}}}; loss = false, delim = "") = loss ? "+" : "-"
decorator(::ElementalScheme{false, <:Groupedisotopomers{<:FormulaChemical}}; loss = false, delim = "") = loss ? "+" : "-"
decorator(::ElementalScheme{false, <:Groupedisotopomers{<:AbstractAdductIon{<:FormulaChemical}}}; loss = false, delim = "") = loss ? "+" : "-"

decorator(::ElementalScheme{true}; loss = false, delim = "") = string(delim, loss ? "Loss_" : "Gain_")
decorator(::ElementalScheme{true, <:FormulaChemical}; loss = false, delim = "") = loss ? "-" : "+"
decorator(::ElementalScheme{true, <:Isotopomers{<:FormulaChemical}}; loss = false, delim = "") = loss ? "-" : "+"
decorator(::ElementalScheme{true, <:Isotopomers{<:AbstractAdductIon{<:FormulaChemical}}}; loss = false, delim = "") = loss ? "-" : "+"
decorator(::ElementalScheme{true, <:Groupedisotopomers{<:FormulaChemical}}; loss = false, delim = "") = loss ? "-" : "+"
decorator(::ElementalScheme{true, <:Groupedisotopomers{<:AbstractAdductIon{<:FormulaChemical}}}; loss = false, delim = "") = loss ? "-" : "+"

charge_repr(c) = c == 0 ? "" : string(abs(c) > 1 ? abs(c) : "", c > 0 ? "+" : "-") 
isotope_repr(isotopes) = isempty(unique_elements(isotopes)) ? "" : string("[", replace(chemicalformula(isotopes; delim = ","), "[" => "", "]" => ""), "]")

chemicalname(cc::Chemical; n = 1, kwargs...) = n > 1 ? string(n, cc.name) : cc.name
function chemicalname(cc::FormulaChemical; n = 1, bracket = false, kwargs...) 
    nm = bracket ? string("[", chemicalformula(cc; kwargs...), "]", charge_repr(charge(cc))) : chemicalformula(cc; kwargs...)
    n > 1 ? string(n, nm) : nm 
end
chemicalname(isobars::Isobars; verbose = true, kwargs...) = (length(chemicalspecies(isobars)) == 1 || verbose) ? string("Isobars[", join(chemicalname.(chemicalspecies(isobars); kwargs...), ", "), "]") : string("Isobars[", chemicalname(first(chemicalspecies(isobars); kwargs...)), ", …]")
chemicalname(isotopomers::Isotopomers; n = 1, kwargs...) = string(chemicalname(chemicalparent(isotopomers); n, kwargs...), isotope_repr(isotopomers.isotopes))
chemicalname(isotopomers::Isotopomers{<:FormulaChemical}; n = 1, kwargs...) = string(chemicalname(chemicalparent(isotopomers); n, bracket = true, kwargs...), isotope_repr(isotopomers.isotopes))

chemicalname(isotopomers::Groupedisotopomers; n = 1, kwargs...) = string(chemicalname(chemicalparent(isotopomers); n, kwargs...), isotopomers.state > 0 ? string("(+", isotopomers.state, ")") : isotopomers.state < 0 ? string("(-", abs(isotopomers.state), ")") : "") 
chemicalname(ct::ChemicalTransition; kwargs...) = join(chemicalname.(chemicaltransition(ct); kwargs...), " -> ")

chemicalname(sch::StructuralElementalScheme; n = 1, kwargs...) = chemicalname(elementalscheme(sch); n, kwargs...)
chemicalname(sch::ElementalScheme{true}; n = 1, loss = false, delim = "", kwargs...) = string(decorator(sch; loss, delim), chemicalname(sch.chemical; n, loss, kwargs...), post_decorator(sch; delim))
chemicalname(sch::ElementalScheme{false}; n = 1, loss = false, delim = "", kwargs...) = string(decorator(sch; loss, delim), chemicalname(sch.chemical; n, loss = !loss, kwargs...), post_decorator(sch; delim))
chemicalname(sch::IsotopomerizedSchema; n = 1, loss = false, bracket = true, delim = "|", kwargs...) = string(chemicalname(chemicalparent(sch); n, loss, bracket, delim, kwargs...), isotope_repr(sch.isotopes))
function chemicalname(sch::ChemicalSchema; n = 1, loss = false, bracket = false, delim = "|", kwargs...) 
    v = [chemicalname(k; n = n * v, loss, bracket = true, delim, kwargs...) for (k, v) in zip(sch.schema, sch.number)]
    predelim = true 
    for (i, s) in enumerate(v)
        if predelim && startswith(s, delim)
            v[i] = s[begin + 1:end]
        end
        predelim = endswith(s, delim)
    end
    s = join(v, "")
    if endswith(s, delim)
        s = s[begin:end - 1]
    end
    bracket ? string("[", s, "]") : s 
end
chemicalname(sch::Groupedisotopomerizedschema; n = 1, loss = false, bracket = true, delim = "|", kwargs...) = string(chemicalname(chemicalparent(sch); n, loss, bracket, delim, kwargs...), sch.state > 0 ? string("(+", sch.state, ")") : sch.state < 0 ? string("(-", abs(sch.state), ")") : "") 

chemicalabbr(isobars::Isobars; verbose = true, kwargs...) = (length(chemicalspecies(isobars)) == 1 || verbose) ? string("Isobars[", join(chemicalabbr.(chemicalspecies(isobars); kwargs...), ", "), "]") : string("Isobars[", chemicalabbr(first(chemicalspecies(isobars); kwargs...)), ", …]")
chemicalabbr(isotopomers::Isotopomers; n = 1, kwargs...) = string(chemicalabbr(chemicalparent(isotopomers); n, kwargs...), isotope_repr(isotopomers.isotopes))
chemicalabbr(isotopomers::Groupedisotopomers; n = 1, kwargs...) = string(chemicalabbr(chemicalparent(isotopomers); n, kwargs...), isotopomers.state > 0 ? string("(+", isotopomers.state, ")") : isotopomers.state < 0 ? string("(-", abs(isotopomers.state), ")") : "") 
chemicalabbr(ct::ChemicalTransition; kwargs...) = join(chemicalabbr.(chemicaltransition(ct); kwargs...), " -> ")

chemicalabbr(sch::StructuralElementalScheme; n = 1, bracket = true, kwargs...) = chemicalabbr(elementalscheme(sch); n, bracket, kwargs...)
chemicalabbr(sch::ElementalScheme{true, <:Electron}; n = 1, bracket = true, loss = false, kwargs...) = ""
chemicalabbr(sch::ElementalScheme{false, <:Electron}; n = 1, bracket = true, loss = false, kwargs...) = ""
chemicalabbr(sch::ElementalScheme{true}; n = 1, bracket = true, loss = false, kwargs...) = string(loss ? "-" : "+", chemicalabbr(sch.chemical; n, bracket, loss, kwargs...))
chemicalabbr(sch::ElementalScheme{false}; n = 1, bracket = true, loss = false, kwargs...) = string(loss ? "+" : "-", chemicalabbr(sch.chemical; n, bracket, loss = !loss, kwargs...))
chemicalabbr(sch::IsotopomerizedSchema; loss = false, bracket = true, n = 1, kwargs...) = string(chemicalabbr(chemicalparent(sch); n, loss, bracket, kwargs...), isotope_repr(sch.isotopes))
function chemicalabbr(sch::ChemicalSchema; loss = false, bracket = true, n = 1, kwargs...) 
    s = join([chemicalabbr(k; n = n * v, loss, bracket = false, kwargs...) for (k, v) in zip(sch.schema, sch.number)], "")
    bracket ? string("[", s, "]", charge_repr(charge(sch; loss))) : s 
end
chemicalabbr(sch::Groupedisotopomerizedschema; n = 1, loss = false, bracket = true, kwargs...) = string(chemicalname(chemicalparent(sch); n, loss, bracket, kwargs...), sch.state > 0 ? string("(+", sch.state, ")") : sch.state < 0 ? string("(-", abs(sch.state), ")") : "") 

chemicalsmiles(isobars::Isobars; kwargs...) = chemicalsmiles(chemicalentity(isobars); kwargs...)
chemicalsmiles(isotopomers::Isotopomers; kwargs...) = chemicalsmiles(chemicalparent(isotopomers); kwargs...)
chemicalsmiles(isotopomers::Groupedisotopomers; kwargs...) = chemicalsmiles(chemicalentity(isotopomers); kwargs...)
chemicalsmiles(ct::ChemicalTransition; kwargs...) = chemicalsmiles(chemicalentity(ct); kwargs...)

function Base.show(io::IO, cc::AbstractChemical)
    nm = chemicalname(cc)
    print(io, isempty(nm) ? chemicalformula(cc) : nm)
end

function Base.show(io::IO, sch::AbstractScheme) 
    print(io, chemicalname(sch))
end

repr_ri(ri::IntervalSet) = isempty(ri) ? "∅" : join([repr_ri(i) for i in ri.items], "∪")
function repr_ri(ri::Interval{T, L, R}) where {T, L, R}
    ff = L == Closed ? "[" : "("
    ll = R == Closed ? "]" : ")"
    f = L == Unbounded ? "-∞" : ri.first
    l = R == Unbounded ? "∞" : ri.last
    string(ff, f, ", ", l, ll)
end

function Base.show(io::IO, ri::IntervalSet)
    print(io, repr_ri(ri))
end

function Base.show(io::IO, c::Criteria{A, B}) where {A<:IntervalSet, B<:IntervalSet}
    print(io, "Criteria{IntervalSet, IntervalSet}(")
    print(io, repr_ri(c.aval), ", ")
    print(io, repr_ri(c.rval), ")")
end

function Base.show(io::IO, c::Criteria{A, B}) where {A<:Missing, B<:IntervalSet}
    print(io, "Criteria{Missing, IntervalSet}(")
    print(io, "missing, ")
    print(io, repr_ri(c.rval), ")")
end

function Base.show(io::IO, c::Criteria{A, B}) where {A<:IntervalSet, B<:Missing}
    print(io, "Criteria{IntervalSet, Missing}(")
    print(io, repr_ri(c.aval), ", ")
    print(io, "missing)")
end

function Base.show(io::IO, spec::Spectrum)
    print(io, "Spectrum of ", round(spec.initial_mass, digits = 4), " ~ ", round(spec.initial_mass + (length(spec.spectrum) - 1) * spec.binsize, digits = 4))
end

function Base.show(io::IO, ci::CoelutingIsobars)
    print(io, "Co-eluting isobars: ", join([measure_name(first(x)) for x in ci.elution], " -> "), " -> ", join([msanalyzer_name(first(x)) for x in ci.msanalyzer], " -> "))
end