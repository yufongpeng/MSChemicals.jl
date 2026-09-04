include(joinpath("isotopologues", "utils.jl"))
include(joinpath("isotopologues", "ms1.jl"))
include(joinpath("isotopologues", "ms2.jl"))
include(joinpath("isotopologues", "msn.jl"))

"""
    Isotopologues(chemical; abundance = 1, abtype = :max, threshold = rcrit(1e-4), precise = false, sort = false)
    Isotopologues(formula_name; chemicalparser, abundance = 1, abtype = :max, threshold = rcrit(1e-4), precise = false, sort = false)
    Isotopologues(chemicaltransition; abundance = 1, abtype = :max, threshold = rcrit(1e-4), precise = false, sort = false)
    Isotopologues(formula_name_pair; chemicalparser, abundance = 1, abtype = :max, threshold = rcrit(1e-4), precise = false, sort = false)
    Isotopologues(tbl; threading = nothing, chemicalparser, threshold = rcrit(1e-4), kwargs...)
    Isotopologues(chemicals; kwargs...)

A `Table` of isotopologues of 
* `chemical::AbstractChemical`: a single chemical entity.
* `formula_name::AbstractString`: a chemical formula or name.
* `chemicaltransition::ChemicalTransition`: MSⁿ transition.
* `formula_name_pair::Pair`: MSⁿ transition of formulas or names.
* `tbl::Table`: multiple chemicals in column `Chemical` with abundance in column `Abundance1`, `Abundance2`, ... (optional). 
* `chemicals::Vector`: multiple chemicals.

This function is similar to [`TandemIsotopologues`](@ref); the key difference is that it is iterative and abundance is calculated in the last stage. It performs faster for multiple MS stages and abundance is normalized and filtered at the end.
Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. 

For MSⁿ transition transition, product can be any scheme, including `<:AbstractStructuralScheme`, `<:AbstractStructuralScheme` or formula starting with `-` or `+` (See [`ChemicalExpressionParser`](@ref) for valid string). 

# Keyword Arguments
* `chemicalparser::AbstractChemicalParser`: parser for `formula_name` or `formula_name_pair`. The default parser is `ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0))` (See [`ChemicalTransitionParser`](@ref)).
* `abundance` sets the abundance of the isotope specified by `abtype`. When the input is MS/MS transition, this sets the abundanc of detected chemical. It can also be column `Abundance` of `tbl`.
* `abtype`.
    * `:max`: the most abundant isotopologue.
    * `:input`: the input isotopologue.
    * `:list`: sum of listed isotopologues.
    * `:total`: sum of total isotopologues.
* `threshold` can be a number or [`Criteria`](@ref), representing the lower limit of abundance (absolute and/or relative to maximal value of each spectrum). 
* `threading`: force to use multiple threads (`true`) or single thread (`false`); `nothing` lets the program determine. 
* `precise`: whether using `Bigfloat` for abundance computation.
* `sort`: whether to sort the results by mass.

!!! note "Special precaution for applying to MSⁿ transition"
    Product must come from a single part or mutiple non-overlapping parts of precursor. Isobaric or isomeric products are not considered. For instance, 
    * PC 18:0/18:0 and fatty acyl 18:0 fragment is valid because two fatty acids are independent and identical. 
    * PC 18:0[D5]/18:0 and fatty acyl 18:0[D5] fragment is valid but the contribution of another fatty acid 18:0 is not considered and addional computation of this pair and summation with knowledge of fragmentation efficiency are required for the correct abundances. 
    * PC 18:1/18:0 and fatty acyl 18:0 fragment is also valid but requires additional computation of isobaric contribution of another fatty acid 18:1. 

!!! note "Special precaution for applying to MSⁿ transition with chemical gain"
    Any intermediate [`ChemicalGain`](@ref) are not allowed.
"""
Isotopologues(input_chemical::AbstractChemical; 
        chemicalparser = ChemicalTransitionParser(),
        id = (1, ), 
        abundance = 1, 
        abtype = Max(), 
        threshold = rcrit(1e-4),
        precise = false,
        iter = false,
        sort = false
    ) = 
    Table(Isotopologues_ms1(Val(precise), input_chemical; chemicalparser, id, abundance, abtype, threshold, iter, sort))

function Isotopologues_ms1(precise::Val, input_chemical::AbstractChemical; 
        chemicalparser = ChemicalTransitionParser(),
        id = (1, ), 
        abundance = 1, 
        abtype = Max(), 
        threshold = rcrit(1e-4),
        iter = false,
        sort = false
    ) 
    it = isotopologues_elements_ms1(precise, input_chemical, chemicalelements(input_chemical), first(abundance), abtype, threshold, iter, sort)
    net_charge = charge(input_chemical)
    if iter && net_charge == 0 
        (; ID = [id for _ in eachindex(it.Chemical)], Chemical = it.Chemical, Mmi1 = it.Mass, Abundance1 = it.Abundance, Preab = it.Preab)
    elseif iter
        (; ID = [id for _ in eachindex(it.Chemical)], Chemical = it.Chemical, MZ1 = it.Mass, Abundance1 = it.Abundance, Preab = it.Preab)
    elseif net_charge == 0 
        (; ID = [id for _ in eachindex(it.Chemical)], Chemical = it.Chemical, Mmi1 = it.Mass, Abundance1 = it.Abundance)
    else
        (; ID = [id for _ in eachindex(it.Chemical)], Chemical = it.Chemical, MZ1 = it.Mass, Abundance1 = it.Abundance)
    end
end

Isotopologues(input_chemical::AbstractString; 
        chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
        id = (1, ), 
        abundance = 1, 
        abtype = Max(), 
        threshold = rcrit(1e-4),
        precise = false,
        iter = false,
        sort = false) = 
    Isotopologues(parse_chemical(chemicalparser, input_chemical); id, abundance, abtype, threshold, precise, iter, sort)

function Isotopologues(ct::ChemicalTransition; 
        chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
        id = ntuple(i -> 1, msstage(ct)), 
        abundance = 1, 
        abtype = Max(), 
        threshold = rcrit(1e-4),
        precise = false,
        iter = false,
        sort = false
    ) 
    abtype = abtyped(abtype)
    msstage(ct) == 1 && return Isotopologues(analyzedchemical(ct); id, abundance, abtype, threshold, precise, iter, sort)
    abundance = float(first(abundance))
    trans = chemicaltransition(ct)
    for c in @view trans[begin:end - 1]
        isgainscheme(c) && throw(ArgumentError("Chemical gain can only be in the last stage for `Isotopologues`; use `TandemIsotopologues` instead."))
    end
    precursor_info = serieschemicaldata(ct)
    if isgainscheme(last(trans))
        if length(trans) > 2
            @inbounds for (pre, post) in @views zip(precursor_info[begin:end - 2], precursor_info[begin + 1:end - 1])
                loss_elements!(last(pre), last(post))
            end
        end
        precursor_info = vcat(precursor_info[begin:end - 1], (last(trans), last(trans), dictionary_elements(chemicalelements(last(trans); loss = false))))
    else
        @inbounds for (pre, post) in @views zip(precursor_info[begin:end - 1], precursor_info[begin + 1:end])
            loss_elements!(last(pre), last(post))
        end
    end
    element_vp_vec = Vector{Vector}(undef, length(precursor_info))
    msfix_vec = Vector{float(Int)}(undef, length(precursor_info))
    @inbounds for (i, info) in enumerate(precursor_info)
        element_vp_vec[i], msfix_vec[i] = get_element_fixmass(unique_elements(last(info)))
    end
    precursor = first.(precursor_info)
    it = isotopologues_elements_msn(Val(precise), element_vp_vec, msfix_vec, abundance, abtype, threshold, sort) 
    if isgainscheme(last(trans)) 
        @inbounds for m in it.Mass
            gain = m[end]
            m[end] = m[end - 1] 
            for i in eachindex(m[begin:end - 1])
                m[i] = m[i] - gain
            end
        end
        universal_loss_elements!(it.Element)
    end
    chemical = seriesisotopomerize(precursor, it.Element)
    net_charge = charge.(precursor)
    abs_charge = [max(1, abs(x)) for x in net_charge]
    colab = Symbol(string("Abundance", length(precursor_info))) 
    colmz = all(==(0), net_charge) ? [Symbol(string("Mmi", i)) for i in eachindex(precursor_info)] : [Symbol(string("MZ", i)) for i in eachindex(precursor_info)]
    Table(; 
        ID = [id for _ in eachindex(chemical)], 
        Chemical = chemical, 
        (colmz[i] => net_charge[i] == 0 ? [m[i] for m in it.Mass] : [m[i] / abs_charge[i] + (net_charge[i] < 0) * ME for m in it.Mass] for i in eachindex(precursor_info))..., 
        (colab => it.Abundance, )...
    ) 
end

function Isotopologues(input_chemical::Pair; 
        chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
        id = nothing, 
        abundance = 1, 
        abtype = Max(), 
        threshold = rcrit(1e-4),
        precise = false,
        iter = false,
        sort = false
    ) 
    ct = parse_chemical(chemicalparser, input_chemical)
    Isotopologues(ct; id = isnothing(id) ? ntuple(i -> 1, msstage(ct)) : id, abundance, abtype, threshold, precise, iter, sort)
end

function Isotopologues(mztable::Table; threading = nothing, chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)), threshold = rcrit(1e-4), kwargs...)
    :Chemical in propertynames(mztable) || throw(ArgumentError("No column `Chemical` in input table."))
    mztable = Table(mztable; Chemical = parse_chemical.(Ref(chemicalparser), mztable.Chemical; kwargs...))
    allequal(msstage, mztable.Chemical) || throw(ArgumentError("Chemicals have to be in the same MS stage."))
    kwargs = Dict(kwargs...)
    vec_key = Symbol[]
    colab = allcolnum(propertynames(mztable), "Abundance"; error = false)
    if !isempty(colab)
        mztable = Table(mztable; abundance = collect.(getproperties(mztable, colab)))
        push!(vec_key, :abundance)
        delete!(kwargs, :abundance)
        ab = mean(mean, mztable.abundance)
    elseif haskey(kwargs, :abundance)
        mztable = Table(mztable; abundance = vectorize(get(kwargs, :abundance, nothing), length(mztable)))
        push!(vec_key, :abundance)
        delete!(kwargs, :abundance)
        ab = mean(mean, mztable.abundance)
    else
        ab = 1.0
    end
    if !in(:ID, propertynames(mztable))
        msn = msstage(first(mztable.Chemical)) - 1
        mztable = Table(mztable; ID = [(i, ntuple(j -> 1, msn)...) for i in eachindex(mztable)])
    end
    rn = min(length(mztable), Threads.nthreads())
    if isnothing(threading)
        s = mean(length(chemicalelements(x)) for x in mztable.Chemical)
        b = mean(mean(last, chemicalelements(x)) for x in mztable.Chemical)
        b = sum(b ^ (1/x) for x in 1:msstage(first(mztable.Chemical)))
        # println((0.02b + 0.2sqrt(-2 * b * log2(min(1, minimum(makecrit_value(crit(threshold), ab)) / ab)))) ^ 1.2s * (rn - 1))
        threading = (0.02b + 0.2sqrt(-2b * log2(min(1, minimum(makecrit_value(crit(threshold), ab)) / ab)))) ^ 1.2s * (rn - 1) > 1e6
    end
    if threading
        t = Vector{Table}(undef, length(mztable))
        Threads.@threads for i in eachindex(t)
            t[i] = Isotopologues(mztable.Chemical[i]; kwargs..., id = mztable.ID[i], [k => getproperty(mztable, k)[i] for k in vec_key]..., threshold)
        end
        tbl = Table(; (p => ChainedVector(getproperty.(t, p)) for p in propertynames(t[1]))...)
    else
        t = [Isotopologues(r.Chemical; kwargs..., id = r.ID, [k => getproperty(r, k) for k in vec_key]..., threshold) for r in mztable]
        tbl = Table(; (p => ChainedVector(getproperty.(t, p)) for p in propertynames(t[1]))...)
    end
    # spectrum specific threshold ?
    colab = lastcolnum(propertynames(tbl), "Abundance"; error = false)
    ab = getproperty(tbl, colab)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>=(abundance_cutoff), ab)
    tbl[id]
end

Isotopologues(v::Vector; kwargs...) = Isotopologues(Table(; Chemical = v); kwargs...) 
Isotopologues(::Isobars; kwargs...) = throw(ArgumentError("`Isobars` is not supported by `Isotopologues.`"))
Isotopologues(::Isotopomers; kwargs...) = throw(ArgumentError("`Isotopomers` is not supported by `Isotopologues.`"))

"""
    TandemIsotopologues(chemical; chemicalparser, product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4), precise = false, sort = false)
    TandemIsotopologues(formula_name; chemicalparser, product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4), precise = false, sort = false)
    TandemIsotopologues(chemicaltransition; chemicalparser, product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4), precise = false, sort = false)
    TandemIsotopologues(formula_name_pair; chemicalparser, product = nothing, proportion = nothing, abundance = 1, abtype = :max, threshold = rcrit(1e-4), precise = false, sort = false)
    TandemIsotopologues(tbl; threading = nothing, chemicalparser, threshold = rcrit(1e-4), kwargs...)
    TandemIsotopologues(chemicals; kwargs...)

A `Table` of isotopologues of the following chemicals and their products. 
* `chemical::AbstractChemical`: a single chemical entity.
* `formula_name::AbstractString`: a chemical formula or name.
* `chemicaltransition::ChemicalTransition`: MSⁿ transition.
* `formula_name_pair::Pair`: MSⁿ transition of formulas or names. 
* `tbl::Table`: multiple chemicals in column `Chemical` with abundance in column `Abundance1`, `Abundance2`, ... (optional). 
* `chemicals::Vector`: multiple chemicals.

This function is similar to `Isotopologues`; the key difference is that it is recursive and abundance is calculated from the beginning. It generally performs slightly slower for multiple MS stages and abundance is normalized in the first stage and filtered in all stages.
Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. 

For MSⁿ transition transition, product can be any scheme, including `<:AbstractStructuralScheme`, `<:AbstractStructuralScheme` or formula starting with `-` or `+` (See [`ChemicalExpressionParser`](@ref) for valid string). 

# Keyword Arguments
* `chemicalparser::AbstractChemicalParser`: parser for `formula_name`, `formula_name_pair` or `product`. The default parser is `ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0))` (See [`ChemicalTransitionParser`](@ref)).
* `abundance` sets the abundance of the precursor isotope specified by `abtype`. It can be a vector when the input is MS/MS transition, and abundances are matched to MS stages from the end (the last element matches to the last MS stage). 
If the length of abundance is smaller, the remaining elemets are filled using `transmission`. It can also be column `Abundance` of `tbl`.
* `abtype`.
    * `:max`: the most abundant isotopologue.
    * `:input`: the input isotopologue.
    * `:list`: sum of listed isotopologues.
    * `:total`: sum of total isotopologues.
* `threshold` can be a number or [`Criteria`](@ref), representing the lower limit of abundance (absolute and/or relative to maximal value of each spectrum). 
* `product::Vector`: product chemicals. It can also be column `Product` of `tbl`.
* `transmission`: transmission rate between precursors (MS/MS transition). It is utlized when all elements of `abundance` are used out. It can also be column `Transmission` of `tbl`.
* `proportion::Vector`: proportion of fragmentation relative to precursor. It can also be column `Proportion` of `tbl`.
* `threading`: force to use multiple threads (`true`) or single thread (`false`); `nothing` lets the program determine. 
* `precise`: whether using `Bigfloat` for abundance computation.
* `sort`: whether to sort the results by mass.

!!! note "Setting abundance ≠ particular isotopologue abundance"
    Setting abundance does not guarantee the equality of abundance of a particular isotopologue in each MS stage, but the abundance sum of isotopologues derived from the particular input precursor isotopologue. 

!!! note "Special precaution for applying to MSⁿ transition"
    Product must come from a single part or mutiple non-overlapping parts of precursor. Isobaric or isomeric products are not considered. For instance, 
    * PC 18:0/18:0 and fatty acyl 18:0 fragment is valid because two fatty acids are independent and identical. 
    * PC 18:0[D5]/18:0 and fatty acyl 18:0[D5] fragment is valid but the contribution of another fatty acid 18:0 is not considered and addional computation of this pair and summation with knowledge of fragmentation efficiency are required for the correct abundances. 
    * PC 18:1/18:0 and fatty acyl 18:0 fragment is also valid but requires additional computation of isobaric contribution of another fatty acid 18:1. 

!!! note "Special precaution for applying to MSⁿ transition with chemical gain"
    After any [`ChemicalGain`](@ref), the subsequent products are considered randomly fragmented from the gained precursor without considering any structure introduced by chemical gain.
"""
function TandemIsotopologues(input_chemical::AbstractChemical; 
            chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
            abundance = 1, 
            transmission = 1, 
            abtype = Max(), 
            threshold = rcrit(1e-4), 
            id = nothing, 
            precursor_table = nothing, 
            precursor_info = nothing, 
            product = nothing, 
            product_info = nothing, 
            proportion = nothing,
            precise = false,
            iter = false,
            sort = false
        ) 
    id = if isnothing(id) 
        isnothing(precursor_table) ? ntuple(x -> 1, msstage(input_chemical)) : first(precursor_table.ID)
    else
        id 
    end
    end_stage = isnothing(product) || isempty(product)
    end_stage && length(id) == 1 && return Isotopologues(input_chemical; chemicalparser, abundance, abtype, threshold, id, precise, iter, sort)
    precursor_info = isnothing(precursor_info) ? serieschemicaldata(input_chemical) : precursor_info
    abundance = vectorize(abundance)
    if length(abundance) > length(precursor_info)
        abundance = abundance[begin:begin + length(precursor_info) - 1]
    elseif length(abundance) < length(precursor_info)
        abundance = vcat(reverse([first(abundance) / transmission ^ i for i in 1:(length(precursor_info) - length(abundance))]), abundance)
    end
    if isnothing(precursor_table)
        precursor_table = TandemIsotopologues_precursor(Val(precise), precursor_info, id, abundance; chemicalparser, abtype, threshold, iter = end_stage ? iter : true, sort)
    end
    end_stage && return Table(precursor_table; Chemical = ChemicalTransition.(precursor_table.Chemical))
    precursor_sch, precursor, element_precursor = last(precursor_info)
    product = parse_chemical.(Ref(chemicalparser), product)
    all(x -> x isa AbstractScheme || msstage(x) < 2, product) || throw(ArgumentError("Products should not be MS/MS pairs."))
    proportion = if isnothing(proportion) 
        [1/length(product) for x in eachindex(product)]
    else
        proportion = vectorize(proportion)
    end
    length(proportion) == length(product) || throw(ArgumentError("The length of `proportion` does not mactch the length of `product`"))
    product_info = if isnothing(product_info) || isempty(product_info)
        map(inputsch -> detectedchemicaldata(precursor, inputsch), product) 
    else
        vectorize(product_info)
    end
    length(product_info) == length(product) || throw(ArgumentError("The length of `product_info` does not mactch the length of `product`"))
    colab = Symbol(string("Abundance", length(id)))
    element_precursor_dictionary = get_element_dictionary(element_precursor)
    isotopes_precursor = map(detectedisotopes, precursor_table.Chemical)
    major_precursor_dictionary = Dict(major_isotope(k) => v for (k, v) in element_precursor_dictionary)
    el = map(isotopes_precursor) do c
        major_minor_precursor_dictionary = copy(major_precursor_dictionary)
        for (k, v) in c
            p = major_isotope(k)
            major_minor_precursor_dictionary[k] = get(major_minor_precursor_dictionary, k, 0) + v
            major_minor_precursor_dictionary[p] = get(major_minor_precursor_dictionary, p, 0) - v
        end
        major_minor_precursor_dictionary
    end
    itp = hasproperty(precursor_table, :Preab) ? (; Element = el, Isotope = isotopes_precursor, Abundance = getproperty(precursor_table, colab), Preab = precursor_table.Preab) : 
        (; Element = el, Isotope = isotopes_precursor, Abundance = getproperty(precursor_table, colab))
    tbls = [TandemIsotopologues_product(Val(precise), precursor_table, itp, element_precursor_dictionary, (id..., i), last(precursor_info), prod_info, prop, last(abundance), Total(), threshold, iter, islossscheme(first(prod_info)), sort; check_product = true) for (i, prop, prod_info) in zip(eachindex(product), proportion, product_info)]
    colab = lastcolnum(propertynames(first(tbls)), "Abundance")
    ab = ChainedVector(getproperty.(tbls, colab))
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = ab .>= abundance_cutoff
    if all(id) 
        tbl = Table(; (p => ChainedVector(getproperty.(tbls, p)) for p in propertynames(first(tbls)))...)
    else
        tbl = Table(; (p => ChainedVector(getproperty.(tbls, p))[id] for p in propertynames(first(tbls)))...)
    end
    Table(tbl; Chemical = [ChemicalTransition(x...) for x in tbl.Chemical])
end

function TandemIsotopologues_precursor(precise::Val, precursor_info, id, abundance; 
            chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
            abtype = Max(), 
            threshold = rcrit(1e-4),
            iter = true,
            sort = false
        )
    precursor_table = Isotopologues_ms1(precise, first(precursor_info)[2]; chemicalparser, id = id[begin:begin], abundance = abundance[begin], abtype, threshold, iter = true, sort)
    length(precursor_info) < 2 && return precursor_table
    ip = 1
    iters = trues(length(precursor_info))
    loss = map(islossscheme ∘ first, precursor_info)
    iters[end - 1] = !isgainscheme(first(last(precursor_info)))
    iters[end] = iter
    while ip < lastindex(precursor_info)
        ip += 1
        colab = Symbol(string("Abundance", ip - 1))
        element_precursor_dictionary = get_element_dictionary(last(precursor_info[ip - 1]))
        major_precursor_dictionary = Dict(major_isotope(k) => v for (k, v) in element_precursor_dictionary)
        el = map(precursor_table.Chemical) do c
            major_minor_precursor_dictionary = copy(major_precursor_dictionary)
            for (k, v) in detectedisotopes(c)
                p = major_isotope(k)
                major_minor_precursor_dictionary[k] = get(major_minor_precursor_dictionary, k, 0) + v
                major_minor_precursor_dictionary[p] = get(major_minor_precursor_dictionary, p, 0) - v
            end
            major_minor_precursor_dictionary
        end
        itp = iters[ip - 1] ? (; Element = el, Abundance = getproperty(precursor_table, colab), Preab = precursor_table.Preab) : 
            (; Element = el, Abundance = getproperty(precursor_table, colab))
            # itp = (; Element = el, Isotope = isotopes_precursor, Abundance = getproperty(precursor_table, colab))
        precursor_table = TandemIsotopologues_product(precise, precursor_table, itp, element_precursor_dictionary, id[begin:ip], precursor_info[ip - 1], precursor_info[ip], abundance[ip] / abundance[ip - 1], abundance[ip], Total(), threshold, iters[ip], loss[ip], sort) 
    end
    precursor_table
end

function TandemIsotopologues_product(precise::Val, precursor_table, itp, element_precursor_dictionary, id, precursor_info, product_info, proportion, abundance, abtype, threshold, iter, chemical_loss, sort; check_product = false)
    precursor_sch, precursor, element_precursor = precursor_info
    product_sch, product, element_product = product_info
    nms = length(id)
    if check_product && !isgainscheme(product_sch) 
        for (k, v) in element_product
            i = findfirst(x -> first(x) == k, element_precursor)
            isnothing(i) && v != 0 && throw(ArgumentError(string("Product ", product, " can only contain elements restricted by precursor ", precursor, ".")))
            last(element_precursor[i]) < v && throw(ArgumentError(string("Product ", product, " can only contain elements restricted by precursor ", precursor, ".")))
        end
    end
    net_charge = charge(product)
    gain = isgainscheme(product_sch) 
    if gain
        element_precursor = chemicalelements(product_sch; loss = false)
    end
    it = isotopologues_elements_ms2(precise, precursor_table.Chemical, product_sch, product, itp, element_precursor_dictionary, element_precursor, element_product, abundance, abtype, proportion, threshold, iter, gain, chemical_loss, sort)
    abpre = map(1:nms - 1) do i 
        s = Symbol(string("Abundance", i))
        s => [getproperty(precursor_table, s)[id] for id in it.ID]
    end
    if net_charge == 0 
        mspre = map(1:nms - 1) do i 
            s = Symbol(string("Mmi", i))
            s => [getproperty(precursor_table, s)[id] for id in it.ID]
        end
        if iter 
            (; ID = [id for _ in eachindex(it.Chemical)], Chemical = it.Chemical, mspre..., [Symbol(string("Mmi", nms)) => it.Mass]..., abpre..., [Symbol(string("Abundance", nms)) => it.Abundance]..., Preab = it.Preab)
        else
            (; ID = [id for _ in eachindex(it.Chemical)], Chemical = it.Chemical, mspre..., [Symbol(string("Mmi", nms)) => it.Mass]..., abpre..., [Symbol(string("Abundance", nms)) => it.Abundance]...) 
        end
    else
        mspre = map(1:nms - 1) do i 
            s = Symbol(string("MZ", i))
            s => [getproperty(precursor_table, s)[id] for id in it.ID]
        end
        if iter 
            (; ID = [id for _ in eachindex(it.Chemical)], Chemical = it.Chemical, mspre..., [Symbol(string("MZ", nms)) => it.Mass]..., abpre..., [Symbol(string("Abundance", nms)) => it.Abundance]..., Preab = it.Preab)
        else
            (; ID = [id for _ in eachindex(it.Chemical)], Chemical = it.Chemical, mspre..., [Symbol(string("MZ", nms)) => it.Mass]..., abpre..., [Symbol(string("Abundance", nms)) => it.Abundance]...) 
        end
    end
end

TandemIsotopologues(input_chemical::Pair; 
            chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
            abundance = 1, 
            transmission = 1, 
            abtype = Max(), 
            threshold = rcrit(1e-4), 
            id = nothing, 
            precursor_table = nothing, 
            precursor_info = nothing, 
            product = nothing, 
            product_info = nothing, 
            proportion = nothing,
            precise = false,
            iter = false,
            sort = false
        ) = TandemIsotopologues(parse_chemical(chemicalparser, input_chemical); chemicalparser, abundance, transmission, abtype, threshold, id, precursor_table, precursor_info, product, product_info, proportion, precise, iter, sort)

TandemIsotopologues(input_chemical::AbstractString; 
            chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)),
            abundance = 1, 
            transmission = 1, 
            abtype = Max(), 
            threshold = rcrit(1e-4), 
            id = nothing, 
            precursor_table = nothing, 
            precursor_info = nothing, 
            product = nothing, 
            product_info = nothing, 
            proportion = nothing,
            precise = false,
            iter = false,
            sort = false
        ) = TandemIsotopologues(parse_chemical(chemicalparser, input_chemical); chemicalparser, abundance, transmission, abtype, threshold, id, precursor_table, precursor_info, product, product_info, proportion, precise, iter, sort)

function TandemIsotopologues(mztable::Table; threading = nothing, chemicalparser = ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 0, gain = 0)), threshold = rcrit(1e-4), kwargs...)
    :Chemical in propertynames(mztable) || throw(ArgumentError("No column `Chemical in input table.`"))
    mztable = Table(mztable; Chemical = parse_chemical.(Ref(chemicalparser), mztable.Chemical; kwargs...))
    allequal(msstage, mztable.Chemical) || throw(ArgumentError("Chemicals have to be in the same MS stage."))
    kwargs = Dict(kwargs...)
    vec_key = Symbol[]
    if :Product in propertynames(mztable)
        mztable = Table(mztable; product = mztable.Product)
        push!(vec_key, :product)
        delete!(kwargs, :product)
        np = mean(length(x) for x in mztable.product)
    elseif haskey(kwargs, :product)
        mztable = Table(mztable; product = vectorize(get(kwargs, :product, nothing), length(mztable)))
        push!(vec_key, :product)
        delete!(kwargs, :product)
        np = mean(length(x) for x in mztable.product)
    else
        np = 1.0
    end
    # all(x -> all(y -> msstage(y) < 2, x), mztable.Product) || throw(ArgumentError("Products should not MS/MS pairs."))
    colab = allcolnum(propertynames(mztable), "Abundance"; error = false)
    if !isempty(colab)
        mztable = Table(mztable; abundance = collect.(getproperties(mztable, colab)))
        push!(vec_key, :abundance)
        delete!(kwargs, :abundance)
        ab = mean(mean, mztable.abundance)
    elseif haskey(kwargs, :abundance)
        mztable = Table(mztable; abundance = vectorize(get(kwargs, :abundance, nothing), length(mztable)))
        push!(vec_key, :abundance)
        delete!(kwargs, :abundance)
        ab = mean(mean, mztable.abundance)
    else
        ab = 1.0
    end
    if :Transmission in propertynames(mztable)
        mztable = Table(mztable; transmission = mztable.Transmission)
        push!(vec_key, :transmission)
        delete!(kwargs, :transmission)
    elseif haskey(kwargs, :transmission)
        mztable = Table(mztable; transmission = vectorize(get(kwargs, :transmission, nothing), length(mztable)))
        push!(vec_key, :transmission)
        delete!(kwargs, :transmission)
    end
    if :Proportion in propertynames(mztable)
        mztable = Table(mztable; proportion = mztable.Proportion)
        push!(vec_key, :proportion)
        delete!(kwargs, :proportion)
    elseif haskey(kwargs, :proportion)
        mztable = Table(mztable; proportion = vectorize(get(kwargs, :proportion, nothing), length(mztable)))
        push!(vec_key, :proportion)
        delete!(kwargs, :proportion)
    end
    if !in(:ID, propertynames(mztable))
        msn = msstage(first(mztable.Chemical)) - 1
        mztable = Table(mztable; ID = [(i, ntuple(j -> 1, msn)...) for i in eachindex(mztable)])
    end
    rn = min(length(mztable), Threads.nthreads())
    if isnothing(threading)
        s = mean(length(chemicalelements(x)) for x in mztable.Chemical)
        b = mean(mean(last, chemicalelements(x)) for x in mztable.Chemical)
        b = sum(b ^ (1/x) for x in 1:msstage(first(mztable.Chemical)))
        # println((0.02b + 0.2sqrt(-2 * b * log2(min(1, minimum(makecrit_value(crit(threshold), ab)) / ab)))) ^ 1.2s * (rn - 1))
        threading = np * (0.02b + 0.2sqrt(-2b * log2(min(1, minimum(makecrit_value(crit(threshold), ab)) / ab)))) ^ 1.2s * (rn - 1) > 1e6
    end
    if threading
        t = Vector{Table}(undef, length(mztable))
        Threads.@threads for i in eachindex(t)
            t[i] = TandemIsotopologues(mztable.Chemical[i]; kwargs..., id = mztable.ID[i], [k => getproperty(mztable, k)[i] for k in vec_key]..., threshold)
        end
        tbl = Table(; (p => ChainedVector(getproperty.(t, p)) for p in propertynames(t[1]))...)
    else
        t = [TandemIsotopologues(r.Chemical; kwargs..., id = r.ID, [k => getproperty(r, k) for k in vec_key]..., threshold) for r in mztable]
        tbl = Table(; (p => ChainedVector(getproperty.(t, p)) for p in propertynames(t[1]))...)
    end
    colab = propertynames(tbl)[findlast(x -> startswith(string(x), "Abundance"), propertynames(tbl))]
    ab = getproperty(tbl, colab)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>=(abundance_cutoff), ab)
    tbl[id]
end

TandemIsotopologues(v::Vector; kwargs...) = TandemIsotopologues(Table(; Chemical = v); kwargs...) 
TandemIsotopologues(::Isobars; kwargs...) = throw(ArgumentError("`Isobars` is not supported by `TandemIsotopologues`."))
TandemIsotopologues(::Isotopomers; kwargs...) = throw(ArgumentError("`Isotopomers` is not supported by `TandemIsotopologues`."))

"""
    group_isotopologues(mztable::Table; isotope = "[13C]")

Group isotopologues by isotopomer state based on `isotope`.
"""
function group_isotopologues(mztable::Table; isotope = "[13C]")
    sp = string.(propertynames(mztable))
    colmz = allcolnum(sp, "MZ")
    colab = allcolnum(sp, "Abundance")
    gf = gf_parent_isotope(isotope)
    transitions = chemicaltransition.(mztable.Chemical)
    gid, gt = unique_group_mz_ab(mztable, transitions, colmz, colab, gf)
    # gt = map(gid) do v 
    #     (; [c => mean(getproperty(mztable, c)[v], weights(getproperty(mztable, d)[v])) for (c, d) in zip(colmz, colab)]..., [c => sum(getproperty(mztable, c)[v]) for c in colab]...)
    # end
    chemcial_parent_state = collect(keys(gid))
    chemical_isotopes = [collect(zip([isotopomersisotopes.(x; loss = false) for x in @view transitions[id]]...)) for id in gid]
    chemical_abundance = [[getproperty(mztable, a)[id] for a in colab] for id in gid]
    chemical = [ChemicalSeries([groupedisotopomerize(p..., isotope, collect(i), a) for (p, i, a) in zip(pa, iso, ab)]) for (pa, iso, ab) in zip(chemcial_parent_state, chemical_isotopes, chemical_abundance)]
    Table(Table(; Chemical = chemical), Table(collect(NamedTuple, gt)))
end

"""
    isotopicabundance(chemical::AbstractChemical, total = 1.0; ignore_isotopes = false, precise = false)
    isotopicabundance(formula::AbstractString, total = 1.0; ignore_isotopes = false, precise = false)
    isotopicabundance(elements::Union{<:Vector, <:Dict}, total = 1.0; ignore_isotopes = false, precise = false)

Compute isotopic abundance of `chemical`, `formula`, vector of element-number pairs or dictionary mapping element to number. `total` is the abundance of all isotopologues.

Parent elements are viewed as major isotopes, and isotopic abundances of all elements are considered in computation. 
To compute isotopic abundance of chemicals with all isotopes labeled intentionally and not following natural distribution, set keyword argument `ignore_isotopes` true, and only parent elements are considered.  
"""
isotopicabundance(cc::AbstractChemical, total = 1.0; ignore_isotopes = false, precise = false) = isotopicabundance(chemicalformula(cc), total; ignore_isotopes, precise)
isotopicabundance(formula::AbstractString, total = 1.0; ignore_isotopes = false, precise = false) = isotopicabundance(chemicalelements(formula), total; ignore_isotopes, precise)
isotopicabundance(elements::Vector{<:Pair}, total = 1.0; ignore_isotopes = false, precise = false) = isotopicabundance(Val(precise), unique_elements(elements), total; ignore_isotopes)
isotopicabundance(elements::Pair, total = 1.0; ignore_isotopes = false, precise = false) = isotopicabundance(Val(precise), elements, total; ignore_isotopes)
isotopicabundance(elements::ElementsVector, total = 1.0; ignore_isotopes = false, precise = false) = isotopicabundance(Val(precise), unique_elements(elements), total; ignore_isotopes)
isotopicabundance(elements::Dict, total = 1.0; ignore_isotopes = false, precise = false) = isotopicabundance(Val(precise), collect(elements), total; ignore_isotopes)
isotopicabundance(precise::Val, elements::Dict, total = 1.0; ignore_isotopes = false) = isotopicabundance(precise, collect(elements), total; ignore_isotopes)
isotopicabundance(precise::Val, elements::ElementsVector, total = 1.0; ignore_isotopes = false) = isotopicabundance(precise, collect(elements), total; ignore_isotopes)
function isotopicabundance(precise::Val, elements::Pair, total = 1.0; ignore_isotopes = false)
    update_isotopicabundance(precise, total, elements)
end
function isotopicabundance(precise::Val, elements::Vector{<:Pair}, total = 1.0; ignore_isotopes = false)
    elements = ignore_isotopes ? filter(iselement ∘ first, elements) : elements
    update_isotopicabundance(precise, total, elements)
end
function update_isotopicabundance(precise::Val, total, elements::Pair)
    total *= precise_exp(precise, get(elements_abundance(), first(elements), one(total)), last(elements))
    return_abundance(precise, total)
end
function update_isotopicabundance(precise::Val, total, elements)
    @inbounds for id in groupfind(parent_element ∘ first, elements)
        abundance = [get(elements_abundance(), first(elements[i]), one(total)) for i in id]
        any(==(one(total)), abundance) && continue
        ns = [last(elements[i]) for i in id]
        total *= safe_multinomial(precise, ns) * prod(precise_exp(precise, y, x) for (x, y) in zip(ns, abundance))
    end
    return_abundance(precise, total)
end