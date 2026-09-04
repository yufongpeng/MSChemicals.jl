const MIN_PROPORTION = 1e-3
min_propotion() = MIN_PROPORTION

"""
    detectedchemicaldata(precursor, product)

Detected chemical and dictionary of elements.
"""
function detectedchemicaldata(precursor, product)
    sch = completescheme(precursor, product)
    det = detectedchemical(precursor, sch)
    completeschemechemical(sch), det, unique_elements(Vector{Pair}, chemicalelements(det))
end

"""
    serieschemicaldata(input_chemical; precursor)

Series detected chemical and dictionaries of elements.
"""
function serieschemicaldata(input_chemical)
    sch = AbstractChemicalsSchema[]
    det = AbstractChemical[]
    precursor = nothing
    for c in chemicaltransition(input_chemical) 
        sc = completescheme(precursor, c)
        push!(sch, completeschemechemical(sc))
        push!(det, detectedchemical(precursor, sc))
        precursor = last(det)
    end
    v = map(eachindex(sch)) do i 
        elements_precursor = chemicalelements(det[i])
        sch[i], det[i], unique_elements(Vector{Pair}, elements_precursor)
    end
    for (s, d, e) in v 
        all(x -> last(x) >=0, e) || throw(ArgumentError("Product can only contain elements restricted by precursor."))
    end
    v
end

"""
    seriesisotopomerize(transitions::Vector{<:AbstractChemicalsSchema}, els::Vector{<:Vector{ElementsVector}})

Serial isotomoperize `transitions` with detected isotpic replacements `els`.
"""
function seriesisotopomerize(transitions::Vector{<:AbstractChemicalsSchema}, els::Vector{<:Vector{ElementsVector}})
    @inbounds map(els) do el
        ChemicalTransition(map(enumerate(transitions)) do (i, trans)
            (islossscheme(trans) || isgainscheme(trans)) && i == firstindex(transitions) && throw(ArgumentError("$(typeof(trans)) cannot be input chemical."))
            islossscheme(trans) ? isotopomerize(trans, loss_elements(el[i], el[i - 1])) : isotopomerize(trans, el[i])
        end
        )
    end
end

"""
    maximal_abundance_elements_composition(precise::Val, elements, prev = 1.0) -> Tuple{Vector{<:AbstractFloat}, Vector{Vector{Int}}}

Isotope composition of maximal abundance for each subisotopologues.
"""
function maximal_abundance_elements_composition(precise::Val, elements)
    max_vps = Vector(undef, length(elements))
    max_abs = Vector{typeof(return_abundance(precise, big(0.0)))}(undef, length(elements))
    for (iem, (e, m)) in enumerate(elements)
        m += 1
        xs = elements_isotopes()[e]
        ns = [floor(Int, m * elements_abundance()[x]) for x in xs]
        if sum(ns) >= m
            d = sum(ns) - m + 1
            i = lastindex(ns)
            while d > 0
                if ns[i] > 0
                    ns[i] -= 1
                    d -= 1
                else
                    i -= 1
                end
            end
        elseif sum(ns) < m - 1
            d = m - sum(ns) - 1
            while d > 0
                _, i = findmax([elements_abundance()[x] / (ns[i] + 1) for (i, x) in enumerate(xs)])
                ns[i] += 1 
                d -= 1
            end
        end
        max_vps[iem] = ns
        max_abs[iem] = return_abundance(precise, safe_multinomial(precise, ns) * prod(precise_exp(precise, elements_abundance()[y], x) for (x, y) in zip(ns, xs)))
    end
    max_abs, max_vps
end

"""
    maximal_abundance_elements_composition_check(precise::Val, elements::Vector, abtype) -> Tuple{Vector{<:AbstractFloat}, Vector{Vector{Int}}}

Isotope composition of maximal abundance for each subisotopologues with additional floating point error check.
"""
function maximal_abundance_elements_composition_check(precise::Val, element_vp, abtype)
    first_proportion = isotopicabundance(precise, element_vp)
    max_proportion, max_vp = maximal_abundance_elements_composition(precise, element_vp)
    # max_proportion = isotopicabundance(precise, max_dictionary)
    mp = prod(max_proportion)
    if mp < first_proportion && !isapprox(mp, first_proportion) 
        if precise == Val(false)
            first_proportion = return_abundance(precise, isotopicabundance(Val(true), element_vp))
            max_proportion, max_vp = maximal_abundance_elements_composition(Val(true), element_vp)
            max_proportion = [return_abundance(precise, x) for x in max_proportion]
        elseif mp < first_proportion 
            throw(ArgumentError("The input chemical is too large to estimete isotopic abundance correctly; try set `precise` true."))
        end
    end
    if abtype == Input() && first_proportion / prod(max_proportion) < min_propotion() 
        throw(ArgumentError("Isotopic abundance of input chemical is too small; try use `abtype` other than `Input()`"))
    end
    max_proportion, max_vp
end

"""
    distribute_element(update_fn, precise::Val, prev, e::String, element_precursor::Dict, pn::Int) -> AbstractFloat

For parent element `e`, distribute isotopes from `element_precursor` into product containing `pn` of `e`; then update `prev` using `update_fn`.
"""
function distribute_element(update_fn, precise, prev, e, element_precursor, pn)
    pre = [get(element_precursor, x, 0) for x in elements_isotopes()[e]]
    # pre[1] = get(element_precursor, e, 0)
    en = sum(pre)
    en == 0 && return (prev, zeros(Int, length(pre)))
    # divrem_pro = [divrem(x * pn, en) for x in pre]
    pro = [div(x * pn, en) for x in pre]
    # pro = first.(divrem_pro)
    # spro = sum(pro)
    delta = pn - sum(pro) 
    if delta > 0
        id = sortperm(pre; by = x -> rem(x * pn, en), rev = true)
        i = 1
        while delta > 0 
            pro[id[i]] += 1 
            delta -= 1
            i += 1 
        end
    end
    update_fn(precise, prev, pre, pro), pro
end

"""
    maximal_proportion_composition(precise::Val, element_precursor::Dict, element_product, element_name_loss, prev = 1.0)

Estimate maximal product isotopologue of `element_product` fragmented from `element_precursor`, and compute the vector of numbers of isotopes and proportion relative to all possible isotopologues. 
"""
function maximal_proportion_composition(precise::Val, element_precursor::Dict, element_product, element_name_loss, prev = 1.0)
    max_vps = Vector{Any}(undef, length(element_product))
    for (i, (e, n)) in enumerate(element_product)
        prev, product = distribute_element(update_maximal_proportion, precise, prev, e, element_precursor, n)
        max_vps[i] = product
    end
    for e in element_name_loss
        prev, _ = distribute_element(update_maximal_proportion, precise, prev, e, element_precursor, 0)
    end
    return_abundance(precise, prev), vcat(max_vps...)
end

update_maximal_proportion(::Val{true}, p, pre, pro) = 
    p * multinomial((big(x) for x in pro)...) / multinomial((big(x) for x in pre)...) * multinomial((big(x) - y for (x, y) in zip(pre, pro))...)

function update_maximal_proportion(::Val{false}, p, pre, pro) 
    if check_overflow_multinomial(pre...)
        p * multinomial((big(x) for x in pro)...) / multinomial((big(x) for x in pre)...) * multinomial((big(x) - y for (x, y) in zip(pre, pro))...)
    else
        p * multinomial(pro...) / multinomial(pre...) * multinomial((x - y for (x, y) in zip(pre, pro))...)
    end
end

"""
    maximal_combination_composition(precise::Val, element_precursor::Dict, element_product, element_name_loss, prev = 1.0)

Estimate maximal product isotopologue of `element_product` fragmented from `element_precursor`, and compute the vector of numbers of isotopes and the number of combinations. 
"""
function maximal_combination_composition(precise::Val, element_precursor::Dict, element_product, element_name_loss, prev = 1.0)
    max_vps = Vector{Any}(undef, length(element_product))
    for (i, (e, n)) in enumerate(element_product)
        prev, product = distribute_element(update_maximal_combination, precise, prev, e, element_precursor, n)
        max_vps[i] = product
    end
    for e in element_name_loss
        prev, _ = distribute_element(update_maximal_combination, precise, prev, e, element_precursor, 0)
    end
    return_abundance(precise, prev), vcat(max_vps...)
end

update_maximal_combination(::Val{true}, p, pre, pro) = 
    p * multinomial((big(x) for x in pro)...) * multinomial((big(x) - y for (x, y) in zip(pre, pro))...)

function update_maximal_combination(::Val{false}, p, pre, pro) 
    res = pre .- pro
    if (sum(pro) > sum(res) ? check_overflow_multinomial(pro) : check_overflow_multinomial(res))
        p * multinomial((big(x) for x in pro)...) * multinomial((big(x) for x in res)...)
    else
        p * safe_multinomial(pro) * safe_multinomial(res)
    end
end

"""
    isotopologue_inverse_combination(precise::Val, enumbers, isotopes)

Compute inverse of combination.
"""
function isotopologue_inverse_combination(precise::Val, numbers, isotopes)
    # ignore_isotopes && return 1.0
    return_abundance(precise, mapfoldl(/, groupfind(parent_element, isotopes); init = 1.0) do v 
        safe_multinomial(precise, numbers[v])
    end)
end

"""
    unique_group_mz_ab(mztable, transitions, colmz, colab, gf = nothing)

Group `transitions` by `gf`, create maps from group to id, and extract m/z values and abundance from `mztable`.
"""
function unique_group_mz_ab(mztable, transitions, colmz, colab, gf = nothing)
    if isnothing(gf)
        gfv = [zeros(length(x)) for x in transitions]
    else
        gfv = map(gf, transitions)
    end
    gids = [groupfind(x -> x[begin:i], gfv) for i in eachindex(first(gfv))]
    ma = map(enumerate(gids)) do (i, gid)
        map(gid) do ids
            chemicals = [y[begin:i] for y in @view transitions[ids]]
            uids = ids[[findfirst(x -> x == t, chemicals) for t in unique(chemicals)]]
            ab = getproperty(mztable, colab[i])[uids]
            mean(getproperty(mztable, colmz[i])[uids], weights(ab)), sum(ab)
        end
    end
    gid = last(gids)
    gt = map(keys(gid)) do k
        ks = [k[begin:i] for i in eachindex(gids)]
        ms = [ma[i][ks[i]] for i in eachindex(gids)]
        (; [c => first(m) for (m, c) in zip(ms, colmz)]..., [c => last(m) for (m, c) in zip(ms, colab)]...)
    end
    gid, gt
end

"""
    gf_parent_isotope(isotope = "[13C]")

Create function for grouping isotopologues by parent chemical and isotopomer state defining by `isotope`.
"""
function gf_parent_isotope(isotope = "[13C]")
    isotope_unit = elements_mass()[isotope] - elements_mass()[elements_parents()[isotope]]
    x -> [chemicalparent(m) => isotopomerstate(m; isotope_unit) for m in x]
end

"""
    get_isotopes(elements)

Get all isotopes of each element in `elements`.
"""
get_isotopes(x) = mapreduce(vcat, x) do e 
    elements_isotopes()[first(e)]
end

"""
    get_nisotopes(elements)

Get number isotopes of each element in `elements`.
"""
get_nisotopes(x) = [length(elements_isotopes()[first(e)]) for e in x]

"""
    get_element_dictionary(input_element)

Element dictionary without isotopes.
"""
function get_element_dictionary(input_element)
    element_dictionary = Dict{String, Int}()
    for (e, n) in input_element
        if iselement(e)
            element_dictionary[e] = get(element_dictionary, e, 0) + n
        end
    end
    element_dictionary
end
get_element_dictionary(input_element::Dict) = filter(iselement ∘ first, input_element)

"""
    get_element(input_element) -> Vector{Pair{String, Int}}

Element vector without isotopes.
"""
function get_element(input_element) 
    [e => n for (e, n) in input_element if iselement(e)]
end

"""
    get_fixmass(input_element) -> AbstractFloat

Mass of isotopes.
"""
function get_fixmass(input_element)
    msfix = elements_mass()[""]
    for (e, n) in input_element
        if !iselement(e)
            msfix += elements_mass()[e] * n
        end
    end
    msfix
end

"""
    get_element_dictionary_fixmass(input_element) -> Tuple{Dict,<:AbstractFloat}

Element dictionary without isotopes, and mass of isotopes
"""
function get_element_dictionary_fixmass(input_element)
    element_dictionary = Dict{String, Int}()
    msfix = elements_mass()[""]
    for (e, n) in input_element
        if iselement(e)
            element_dictionary[e] = get(element_dictionary, e, 0) + n
        else
            msfix += elements_mass()[e] * n
        end
    end
    element_dictionary, msfix
end

"""
    get_element_fixmass(input_element) -> Tuple{Vector{Pair{String, Int}},<:AbstractFloat}

Element vector without isotopes, and mass of isotopes
"""
function get_element_fixmass(input_element) 
    msfix = elements_mass()[""]
    for (e, n) in input_element
        if !iselement(e)
            msfix += elements_mass()[e] * n
        end
    end
    [e => n for (e, n) in input_element if iselement(e)], msfix
end

"""
    element_mass_delta(old_element, new_element) -> AbstractFloat

Mass of `new_element` minus mass of `old_element`.
"""
element_mass_delta(old_element, new_element) = elements_mass()[new_element] - elements_mass()[old_element]

"""
    deltammi(elements::Vector{String}, numbers::Vector{Int}) -> AbstractFloat

Sum of mass of each element in `elements` minus mass of their parent elements.
"""
deltammi(isotopes::Vector{String}, v) = sum((elements_mass()[x] - elements_mass()[parent_element(x)]) * n for (x, n) in zip(isotopes, v))

"""
    nmmi(elements::Vector{String}, numbers::Vector{Int}) -> AbstractFloat

Mass of `elements` and `numbers`.
"""
nmmi(isotopes::Vector{String}, v) = sum(elements_mass()[x] * n for (x, n) in zip(isotopes, v))

"""
    update_abundance(precise::Val, prev_abundance, old_element, new_element, nold, nnew) -> AbstractFloat

Update of `prev_abundance` after one `old_element` replaced by `new_element`. `nold` and `nnew` are number of elements before replacing.
"""
function update_abundance(precise::Val, prev_abundance::T, old_element, new_element, nold, nnew) where T
    x = get(elements_abundance(), old_element, one(T))
    y = get(elements_abundance(), new_element, one(T))
    if (x == one(T) || y == one(T))
        prev_abundance
    else 
        update_abundance1(precise, prev_abundance, x, y, nold, nnew)
    end
end

update_abundance1(::Val{true}, prev_abundance, x, y, nold, nnew) = prev_abundance * (big(nold) / (nnew + 1)) * (big(y) / x)
update_abundance1(::Val{false}, prev_abundance, x, y, nold, nnew) = prev_abundance * (nold / (nnew + 1)) * (y / x)

"""
    update_proportion(precise::Val, prev_proportion, nold, nnew) -> AbstractFloat
    
Update of `prev_proportion` after one element replacing. `nold` and `nnew` are number of elements before replacing.
"""
update_proportion(::Val{true}, prev_proportion, nold, nnew) = prev_proportion * (big(nold) / (nnew + 1)) 
update_proportion(::Val{false}, prev_proportion, nold, nnew) = prev_proportion * (nold / (nnew + 1)) 

"""
    update_inverse_proportion(precise::Val, prev_inverse_proportion, nold, nnew) -> AbstractFloat

Update of `prev_inverse_proportion` after one element replacing. `nold` and `nnew` are number of elements before replacing.
"""
update_inverse_proportion(::Val{true}, prev_proportion, nold, nnew) = prev_proportion * (big(nnew + 1) / nold) 
update_inverse_proportion(::Val{false}, prev_proportion, nold, nnew) = prev_proportion * ((nnew + 1) / nold) 