"""
    parse_chemical([parser::AbstractChemicalParser,] name::AbstractString...; kwargs...) 
    parse_chemical([parser::AbstractChemicalParser,] chemical::AbstractChemicalsSchema...; kwargs...) -> AbstractChemicalsSchema
    parse_chemical([parser::AbstractChemicalParser,] chemicals::AbstractVector; kwargs...) 
    parse_chemical([parser::AbstractChemicalParser,] pair::Pair kwargs...) 

Parse chemical name and construct a chemical object using `parser`. The default parser is `ChemicalTransitionParser()` (See [`ChemicalTransitionParser`](ref)).
"""
parse_chemical(name...; kwargs...) = parse_chemical(ChemicalTransitionParser(), name...; kwargs...)
function parse_chemical(chemicalparser::ChemicalParser, name::AbstractString; kwargs...) 
    property = copy(chemicalparser.property)
    ks = first.(property)
    for (k, v) in kwargs
        i = findfirst(==(k), ks)
        if isnothing(i) && !isnothing(v)
            push!(property, k => v)
            push!(ks, k)
        elseif isnothing(i)
            continue
        else
            property[i] = k => v 
        end
    end
    filter!(x -> !isnothing(last(x)), property)
    i = findfirst(==(:formula), first.(property))
    if isnothing(i)
        i = findfirst(==(:elements), first.(property))
        isnothing(i) && throw(ArgumentError("No formula or elements are given"))
    end 
    fe = last(property[i])
    deleteat!(property, i)
    Chemical(name, chemicalelements(fe), property)
end

function parse_chemical(chemicalparser::FormulaChemicalParser, name::AbstractString; kwargs...) 
    property = copy(chemicalparser.property)
    ks = first.(property)
    for (k, v) in kwargs
        i = findfirst(==(k), ks)
        if isnothing(i) && !isnothing(v)
            push!(property, k => v)
            push!(ks, k)
        elseif isnothing(i)
            continue
        else
            property[i] = k => v 
        end
    end
    filter!(x -> !isnothing(last(x)), property)
    FormulaChemical(chemicalelements(name), property)
end

function parse_adduction(chemicalparser, name, core; chemicalgain = nothing, kwargs...)
    if isnothing(core) || isempty(core) 
        corechemical = nothing 
        corecharge = 0
        ncore = 1
    else
        ncore, core = match(r"(\d*)(.*)", core)
        ncore = isempty(ncore) ? 1 : parse(Int, ncore)
        corechemical = isnothing(chemicalgain) ? nothing : get(scheme_abbr(), core, nothing)
        corechemical = isnothing(corechemical) ? parse_chemical(chemicalparser, core; kwargs...) : corechemical
        corecharge = charge(corechemical) * ncore
    end
    ion, nc = match(r"\[(.*)\](\d*[+-]*)", name)
    if isempty(nc)
        net_charge = 0 
    else
        pos = endswith(nc, "+") ? 1 : -1
        nc = nc[begin:end - 1]
        net_charge = (isempty(nc) ? 1 : parse(Int, nc)) * pos
    end
    net_charge -= corecharge
    sch = get(scheme_name(), string("[", ion, "]", charge_repr(net_charge)), nothing)
    isnothing(sch) || return (corechemical, sch, ncore)
    spn = [split(p, "-") for p in split(ion, "+")]
    vc = AbstractScheme[]
    for ss in spn
        if !isempty(first(ss))
            n, fs = match(r"(\d*)(.*)", first(ss))
            n = isempty(n) ? 1 : parse(Int, n)
            fs = ElementalScheme(true, haskey(scheme_abbr(), fs) ? scheme_abbr()[fs] : parse_chemical(chemicalparser, fs; kwargs...))
            for _ in 1:n 
                push!(vc, fs)
            end
        end
        if length(ss) > 1 
            for i in eachindex(ss)[2:end]
                isempty(ss[i]) && continue
                n, fs = match(r"(\d*)(.*)", ss[i])
                n = isempty(n) ? 1 : parse(Int, n)
                fs = ElementalScheme(false, haskey(scheme_abbr(), fs) ? scheme_abbr()[fs] : parse_chemical(chemicalparser, fs; kwargs...))
                for _ in 1:n 
                    push!(vc, fs)
                end
            end
        end
    end
    ne = net_charge - sum(charge, vc; init = 0)
    if ne == 0 && length(vc) == 0
        sch = nothing
    elseif ne == 0 && length(vc) == 1  
        sch = first(vc)
    elseif ne == 1 && length(vc) == 0
        sch = ElementalScheme(false, Electron())
    elseif ne == -1 && length(vc) == 0
        sch = ElementalScheme(true, Electron())
    elseif ne == 0
        sch = ChemicalSchema(vc...)
    elseif ne > 0
        sch = ChemicalSchema(vc..., (ElementalScheme(false, Electron()) for i in 1:ne)...)
    else
        sch = ChemicalSchema(vc..., (ElementalScheme(true, Electron()) for i in ne:(-1))...)
    end
    (corechemical, sch, ncore)
end

assemble_chemical(chemicalgain, corechemical, sch::Nothing, ncore) = ncore > 1 ? ChemicalSchema([ElementalScheme(chemicalgain, corechemical)], [ncore]) : ElementalScheme(chemicalgain, corechemical)
assemble_chemical(chemicalgain::Nothing, corechemical, sch::Nothing, ncore) = ncore > 1 ? throw(ArgumentError("Number of chemicals > 1 is not allowed without any chemical wrapper.")) : corechemical 
assemble_chemical(chemicalgain, corechemical::Nothing, sch, ncore) = throw(ArgumentError("Cannot loss or gain a scheme."))
assemble_chemical(chemicalgain::Nothing, corechemical::Nothing, sch, ncore) = sch 
assemble_chemical(chemicalgain, corechemical, sch, ncore) = ElementalScheme(chemicalgain, AdductIon(corechemical, completescheme(corechemical, sch), ncore))
assemble_chemical(chemicalgain::Nothing, corechemical, sch, ncore) = AdductIon(corechemical, completescheme(corechemical, sch), ncore)

function parse_chemical(chemicalparser::ChemicalExpressionParser, name::AbstractString; precursorcharge = nothing, both = false, kwargs...)
    precursorcharge = isnothing(precursorcharge) ? chemicalparser.charge : precursorcharge
    # process [...]
    if both || chemicalparser.scheme
        # isnothing(core) || return parse_adduction(chemicalparser.chemicalparser, name, core; kwargs...)
        sch = get(scheme_name(), name, nothing)
        isnothing(sch) || return sch 
    end
    # process +/-
    chemicalgain, bracket = isplusminusbracket(name)
    core = nothing
    sch = nothing
    if (both || chemicalparser.scheme) && !isnothing(chemicalgain) && !chemicalgain && !bracket
        chemicalcharge = chemicalparser.loss
        if precursorcharge == 0 
            chemicalcharge = 0
        else
            while (precursorcharge - chemicalcharge) * precursorcharge <= 0 
                chemicalcharge -= sign(precursorcharge)
            end
        end
        chemicalgain = false
        core = name[begin + 1:end]
        sch = string("[]", charge_repr(chemicalcharge))
        both = true
    elseif (both || chemicalparser.scheme) && !isnothing(chemicalgain) && chemicalgain && !bracket
        chemicalcharge = chemicalparser.gain
        if precursorcharge == 0 
            chemicalcharge = 0
        else
            while (precursorcharge + chemicalcharge) * precursorcharge <= 0 
                chemicalcharge += sign(precursorcharge)
            end
        end
        chemicalgain = true
        core = name[begin + 1:end]
        sch = string("[]", charge_repr(chemicalcharge))
        both = true
    elseif (both || chemicalparser.scheme) && !isnothing(chemicalgain)
        both = true
        name = name[begin + 1:end]
    end
    # process chemical part
    if (both || chemicalparser.entity) && startswith(name, "[") 
        core, sch = match(r"\[([^+-]*)(.*\]\d*[+-]*)", name)
        sch = string("[", sch)
    elseif (both || chemicalparser.entity) && isnothing(chemicalgain)
        chemicalcharge = min(abs(precursorcharge), abs(chemicalparser.charge)) * sign(precursorcharge)
        sch = string("[]", charge_repr(chemicalcharge))
        core = name
    elseif (both || chemicalparser.scheme) && startswith(name, r"\[[\]+-]") 
        core = ""
        sch = name
    end
    isnothing(sch) && throw(ArgumentError("Failed to parse \"$name\"."))
    assemble_chemical(chemicalgain, parse_adduction(chemicalparser.chemicalparser, sch, core; kwargs..., chemicalgain)...)
end

function isplusminusbracket(x)
    nb = 0
    x = replace(x, r"\d*[+-]$" => "")
    chemicalgain = startswith(x, "+") ? true : startswith(x, "-") ? false : nothing
    i = firstindex(x)
    start = false
    while i < lastindex(x) 
        if x[i] == '['
            start = true
            nb += 1 
        elseif x[i] == ']' 
            nb -= 1
        end
        if start && nb == 0
            return (chemicalgain, false)
        end
        i = nextind(x, i)
    end
    start ? (chemicalgain, true) : (chemicalgain, false)
end

parse_chemical(::AbstractChemicalParser, cc::AbstractChemicalsSchema...; kwargs...) = ChemicalSeries(cc...)

parse_chemical(::AbstractChemicalParser, cc::AbstractChemicalsSchema, precursorcharge::Union{Int, Nothing}; kwargs...) = (cc, detectedcharge(cc; precursorcharge))
parse_chemical(::AbstractChemicalParser, cc::ChemicalTransition, precursorcharge::Union{Int, Nothing}; kwargs...) = (cc.transition, detectedcharge(cc))
function parse_chemical(chemicalparser::AbstractChemicalParser, name::AbstractString, precursorcharge::Union{Int, Nothing}; kwargs...) 
    chemical = parse_chemical(chemicalparser, name; precursorcharge, kwargs...)
    chemical, detectedcharge(chemical; precursorcharge)
end

parse_chemical(chemicalparser::ChemicalTransitionParser, name::AbstractString...; kwargs...) = parse_chemical(chemicalparser, collect(name); kwargs...)
parse_chemical(chemicalparser::ChemicalTransitionParser, name::AbstractString; kwargs...) = parse_chemical(chemicalparser, split(name, r"\s*->\s*"); kwargs...)
function parse_chemical(chemicalparser::ChemicalTransitionParser, name::AbstractVector; kwargs...) 
    trans = Any[]
    precursorcharge = nothing
    if length(name) < 2
        return parse_chemical(chemicalparser.chemicalparser, first(name); kwargs..., precursorcharge)
    end
    for nm in name
        c, precursorcharge = parse_chemical(chemicalparser.chemicalparser, nm, precursorcharge; kwargs...)
        push!(trans, c)
    end
    ChemicalSeries(trans...)
end

function parse_chemical(chemicalparser::ChemicalTransitionParser, name::Pair; kwargs...) 
    pre, precursorcharge = parse_chemical(chemicalparser, first(name), nothing; kwargs...)
    ChemicalSeries(first(pre_parse_chemical(chemicalparser, pre, last(name), precursorcharge; kwargs...)))
end

function pre_parse_chemical(chemicalparser::ChemicalTransitionParser, pre, name, precursorcharge; kwargs...)
    pre = vectorize(pre) 
    post, precursorcharge = parse_chemical(chemicalparser.chemicalparser, name, detectedcharge(last(pre); precursorcharge); kwargs...)
    vcat(pre, post), precursorcharge
end

function pre_parse_chemical(chemicalparser::ChemicalTransitionParser, pre, name::Pair, precursorcharge; kwargs...) 
    pre, precursorcharge = pre_parse_chemical(chemicalparser, pre, first(name), precursorcharge; kwargs...)
    pre_parse_chemical(chemicalparser, pre, last(name), precursorcharge; kwargs...)
end

function pre_parse_chemical(chemicalparser::ChemicalTransitionParser, pre, name::ChemicalTransition, precursorcharge; kwargs...) 
    post = Any[]
    for nm in chemicaltransition(name)
        c, precursorcharge = parse_chemical(chemicalparser.chemicalparser, nm, precursorcharge; kwargs...)
        push!(post, c)
    end
    vcat(pre, post...), precursorcharge
end


"""
    parse_adduct([adductparser::AbstractAdductParser,] adduct; kwargs...) 

Parse string into `NamedTuple` or `Tuple` using `adductparser`. 
The default `AdductParser` determines the number of core first (`"[2M+...]"` for 2 cores), and then searches string in `scheme_name()`, or parses the string using `adductparser.chemicalparser`.
"""
parse_adduct(adduct; kwargs...) = parse_adduct(AdductParser(), adduct; kwargs...)
parse_adduct(adductparser::AdductParser, adduct::AbstractScheme; kwargs...) = return_adduct(adductparser, adduct, 1)
function parse_adduct(adductparser::AdductParser, adduct::AbstractString; kwargs...) 
    nm, sch = split(adduct, "M"; limit = 2)
    nm = nm[begin + 1:end]
    ncore = isempty(nm) ? 1 : parse(Int, nm)
    sch = string("[", sch)
    adduct = haskey(scheme_name(), sch) ? scheme_name()[sch] : parse_chemical(adductparser.chemicalparser, sch)
    return_adduct(adductparser, adduct, ncore) 
end

return_adduct(::AdductParser{true}, adduct, ncore) = (adduct, ncore)
return_adduct(::AdductParser{false}, adduct, ncore) = (; adduct, ncore)

parse_adduct(x::Tuple; kwargs...) = x
parse_adduct(::AdductParser{true}, x::Tuple; kwargs...) = x
parse_adduct(x::NamedTuple; kwargs...) = x
parse_adduct(::AdductParser{false}, x::NamedTuple; kwargs...) = x
parse_adduct(::AdductParser{true}, x::NamedTuple; kwargs...) = Tuple(x)

# parse_adduct(adduct::AbstractChemicalScheme) = (adduct, 1)
# parse_adduct(::AbstractAdductParser, adduct::AbstractChemicalScheme) = (adduct, 1)

