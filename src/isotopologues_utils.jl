const MIN_PROPORTION = 1e-3
min_propotion() = MIN_PROPORTION

"""
    detectedchemicaldata(precursor, product)

Detected chemical and dictionary of elements.
"""
function detectedchemicaldata(precursor, product)
    sch = completescheme(precursor, product)
    det = detectedchemical(precursor, sch)
    sch, det, unique_elements(Dict, chemicalelements(det))
end

"""
    serieschemicaldata(raw_product; precursor)

Series detected chemical and dictionaries of elements.
"""
function serieschemicaldata(input_chemical)
    sch = AbstractChemicalsSchema[]
    det = AbstractChemical[]
    precursor = nothing
    for c in chemicaltransition(input_chemical) 
        push!(sch, completescheme(precursor, c))
        push!(det, detectedchemical(precursor, last(sch)))
        precursor = last(det)
    end
    v = map(eachindex(sch)) do i 
        elements_precursor = chemicalelements(det[i])
        sch[i], det[i], unique_elements(Dict, elements_precursor)
    end
    for (s, d, e) in v 
        all(>=(0), values(e)) || throw(ArgumentError("Product can only contain elements restricted by precursor."))
    end
    v
end

"""
    seriesisotopomerize(transitions::Vector{<: AbstractChemicalsSchema}, els::Vector{<: Vector{<: Dict}})

Serial isotomoperize `transitions` with detected isotpic replacements `els`.
"""
function seriesisotopomerize(transitions::Vector{<: AbstractChemicalsSchema}, els::Vector{<: Vector{Vector{Pair{String, Int}}}})
    @inbounds map(els) do el
        ChemicalTransition(map(enumerate(transitions)) do (i, trans)
            (islossscheme(trans) || isgainscheme(trans)) && i == firstindex(transitions) && throw(ArgumentError("$(typeof(trans)) cannot be input chemical."))
            islossscheme(trans) ? isotopomerize(trans, loss_elements(el[i], el[i - 1])) : isotopomerize(trans, el[i])
        end
        )
    end
end

# """
#     maximal_abundance(precise::Val, elements::Dict, prev = 1.0) -> AbstractFloat

# Dictionary of elements of maximal abundance.
# """
# function conditional_maximal_abundance(precise::Val, elements, condition, isotopes, addition, prev)
#     e = parent_element(condition)
#     if addition
#         new_p = prev * precise_ratio(precise, elements_abundance()[condition], elements[condition] + 1)
#         prev_p = new_p * precise_ratio(precise, elements[e], elements_abundance()[e])
#         for x in isotopes 
#             max_p = max(prev_p, new_p * precise_ratio(precise, elements[x], elements_abundance()[x]))
#         end
#     else
#         new_p = prev * precise_ratio(precise, elements[condition], elements_abundance()[condition])
#         prev_p = new_p * precise_ratio(precise, elements_abundance()[e], elements[e] + 1)
#         for x in isotopes 
#             max_p = max(prev_p, new_p * precise_ratio(precise, elements_abundance()[x], elements[x] + 1))
#         end
#     end
#     max_p
# end

# precise_ratio(precise::Val{true}, x, y) = big(x) / y
# precise_ratio(precise::Val{false}, x, y) = x / y

"""
    maximal_elements(elements::Dict) -> Dict

Dictionary of elements of maximal abundance.
"""
function maximal_elements(elements)
    max_dictionary = copy(elements)
    for (e, m) in elements
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
        max_dictionary[e] = first(ns)
        i = 1
        while i < lastindex(ns)
            i += 1
            max_dictionary[xs[i]] = ns[i]
        end
    end
    max_dictionary
end

"""
    maximal_abundance(elements::Dict) 

Dictionary of elements of maximal abundance.
"""
function maximal_abundance(precise::Val, elements, prev = 1.0)
    for (e, m) in elements
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
        prev *= safe_multinomial(precise, ns) * prod(precise_exp(precise, elements_abundance()[y], x) for (x, y) in zip(ns, xs))
    end
    return_abundance(precise, prev)
end

"""
    maximal_abundance_elements(elements::Dict) 

Dictionary of elements of maximal abundance.
"""
function maximal_abundance_elements(precise::Val, elements, prev = 1.0)
    max_dictionary = copy(elements)
    for (e, m) in elements
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
        max_dictionary[e] = first(ns)
        i = 1
        while i < lastindex(ns)
            i += 1
            max_dictionary[xs[i]] = ns[i]
        end
        prev *= safe_multinomial(precise, ns) * prod(precise_exp(precise, elements_abundance()[y], x) for (x, y) in zip(ns, xs))
    end
    return_abundance(precise, prev), max_dictionary
end

return_abundance(::Val{true}, x) = x 
return_abundance(::Val{false}, x) = convert(float(Int), x) 

"""
    distribute_element_update!(update_fn, precise, prev, e, isotopes, element_precursor, element_product, max_dict = nothing)

Distribute isotopes from `element_precursor` into `element_product` for parent element `e` and isotopes `isotopes`; then update max_dict and update `prev` using `update_fn`.
"""
function distribute_element_update!(update_fn, precise, prev, e, isotopes, element_precursor, element_product, max_dict = nothing)
    pre = [get(element_precursor, x, 0) for x in isotopes]
    en = sum(pre)
    en == 0 && return prev
    pn = get(element_product, e, 0)
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
    update_max_dict!(max_dict, e, isotopes, pro)
    update_fn(precise, prev, pre, pro)
end

update_max_dict!(max_dict::Nothing, e, isotopes, pro) = nothing
function update_max_dict!(max_dict::Dict, e, isotopes, pro) 
    for (x, n) in zip(isotopes, pro) 
        # x == e && continue
        max_dict[x] = get(max_dict, x, 0) + n
        max_dict[e] = get(max_dict, e, 0) - n
    end
end

"""
    maximal_proportion(precise::Val, element_precursor::Dict, element_product::Dict, prev = 1.0) -> AbstractFloat

Estimate maximal product isotopologue of `element_product`, and compute the proportion relative to all possible isotopologues fragmented from `element_precursor`. 
"""
function maximal_proportion(precise::Val, element_precursor, element_product, prev = 1.0)
    for (e, v) in pairs(group(parent_element, keys(element_precursor))) 
        prev = distribute_element_update!(update_maximal_proportion, precise, prev, e, v, element_precursor, element_product)
    end
    return_abundance(precise, prev)
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
    maximal_proportion_elements(precise::Val, element_precursor::Dict, element_product::Dict, prev = 1.0)

Estimate maximal product isotopologue of `element_product`, and compute `Dictionary` of elements and proportion relative to all possible isotopologues fragmented from `element_precursor`. 
"""
function maximal_proportion_elements(precise::Val, element_precursor::Dict, element_product::Dict, prev = 1.0)
    first_element_product_dictionary = copy(element_product)
    for (e, v) in pairs(group(parent_element, keys(element_precursor))) 
        prev = distribute_element_update!(update_maximal_proportion, precise, prev, e, v, element_precursor, element_product, first_element_product_dictionary)
    end
    return_abundance(precise, prev), first_element_product_dictionary
end

"""
    maximal_combination(precise::Val, element_precursor::Dict, element_product::Dict, prev = 1.0) -> AbstractFloat

Estimate maximal product isotopologue of `element_product`, and compute the number of combinations. 
"""
function maximal_combination(precise::Val, element_precursor::Dict, element_product::Dict, prev = 1.0)
    for (e, v) in pairs(group(parent_element, keys(element_precursor))) 
        prev = distribute_element_update!(update_maximal_combination, precise, prev, e, v, element_precursor, element_product)
    end
    return_abundance(precise, prev)
end

"""
    maximal_combination_elements(precise::Val, element_precursor::Dict, element_product::Dict, prev = 1.0)

Estimate maximal product isotopologue of `element_product`, and compute `Dictionary` of elements and the number of combinations. 
"""
function maximal_combination_elements(precise::Val, element_precursor, element_product, prev = 1.0)
    first_element_product_dictionary = copy(element_product)
    for (e, v) in pairs(group(parent_element, keys(element_precursor))) 
        prev = distribute_element_update!(update_maximal_combination, precise, prev, e, v, element_precursor, element_product, first_element_product_dictionary)
    end
    return_abundance(precise, prev), first_element_product_dictionary
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
    isotopologue_inverse_combination(precise::Val, elements::Dict)
    isotopologue_inverse_combination(precise::Val, element_precursor::Dict, element_product::Dict, swap::Vector{String})

Compute inverse of combination.
"""
function isotopologue_inverse_combination(precise::Val, elements::Dict)
    # ignore_isotopes && return 1.0
    return_abundance(precise, mapfoldl(/, group(parent_element ∘ first, elements); init = 1.0) do v 
        safe_multinomial(precise, last.(v))
    end)
end

function isotopologue_inverse_combination(precise::Val, element_precursor::Dict, element_product::Dict, swap::Vector{String})
    return_abundance(precise, mapfoldl(/, pairs(group(parent_element ∘ first, element_precursor)); init = 1.0) do (e, v) 
        if e in swap 
            safe_multinomial(precise, (n - get(element_product, i, 0) for (i, n) in v)...)
        else
            safe_multinomial(precise, (get(element_product, i, 0) for (i, n) in v)...)
        end
    end)
end

"""
    loss_inverse_combination(precise::Val, element_precursor::Dict, element_product::Dict, swap::Vector{String})

Compute inverse of combination for chemical loss.
"""
function loss_inverse_combination(precise::Val, element_precursor::Dict, element_product::Dict, swap::Vector{String})
    return_abundance(precise, mapfoldl(/, pairs(group(parent_element ∘ first, element_precursor)); init = 1.0) do (e, v) 
        if e in swap 
            safe_multinomial(precise, (get(element_product, i, 0) for (i, n) in v)...)
        else
            safe_multinomial(precise, (n - get(element_product, i, 0) for (i, n) in v)...)
        end
    end)
end

"""
    unique_group_mz_ab(mztable, transitions, colmz, colab, gf = nothing)

Group `transitions` by `gf`, create maps from group to id, and extract m/z values and abundance from `matable`.
"""
function unique_group_mz_ab(mztable, transitions, colmz, colab, gf = nothing)
    if isnothing(gf)
        gfv = [zeros(length(x)) for x in transitions]
    else
        gfv = map(gf, transitions)
    end
    gids = [groupfind(x -> x[begin:i], gfv) for i in eachindex(first(gfv))]
    ma = map(enumerate(gids)) do (i, gid)
        map(gid) do v 
            c = [y[begin:i] for y in transitions[v]]
            u = v[[findfirst(x -> x == t, c) for t in unique(c)]]
            ab = getproperty(mztable, colab[i])[u]
            mean(getproperty(mztable, colmz[i])[u], weights(ab)), sum(ab)
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
    swap_elements(product::Dict, precursor::Dict, swap::Vector{String})

Swap elements of `products` and elements of residuals of `precursor` and `product` for eleemnts in `swap`.
"""
swap_elements(product::Dict, precursor::Dict, swap::Vector{String}) = isempty(swap) ? get_isotope_vec(product) : [k => nswap(k, v, swap, product) for (k, v) in precursor if !iselement(k)]

nswap(k, v, swap, product) = parent_element(k) in swap ? v - get(product, k, 0) : get(product, k, 0) 

"""
    get_isotope_vec(input_element::Vector{Pair{String, Int}})
    get_isotope_vec(input_element::Dict{String, Int})

Create isotopic replacement vectors.
"""
get_isotope_vec(input_element::Vector{Pair{String, Int}}) = filter(x -> !iselement(first(x)), input_element)
get_isotope_vec(input_element::Dict{String, Int}) = [k => v for (k, v) in input_element if !iselement(k)]

"""
    get_element_dictinonary(input_element)

Element dictionary without isotopes.
"""
function get_element_dictinonary(input_element)
    element_dictionary = Dict{String, Int}()
    for (e, n) in input_element
        if iselement(e)
            element_dictionary[e] = get(element_dictionary, e, 0) + n
        end
    end
    element_dictionary
end
get_element_dictinonary(input_element::Dict) = filter(iselement ∘ first, input_element)

"""
    get_fixmass(input_element)

Mass of isotopes
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
    get_element_dictinonary_fixmass(input_element)

Element dictionary without isotopes, and mass of isotopes
"""
function get_element_dictinonary_fixmass(input_element)
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
    element_isotope_pairs(element_dictionary; sort = true)

Elements and isotopes pairs.
"""
function element_isotope_pairs(element_dictionary; sort = true)
    element_isotope_pair = mapreduce(vcat, collect(keys(element_dictionary))) do e
        v = get(elements_isotopes(), e, nothing)
        isnothing(v) && return Pair{String, String}[]
        length(v) < 1 && return Pair{String, String}[]
        map(v[begin + 1:end]) do x
            (e, x)
        end
    end
    sort ? sort!(element_isotope_pair; by = pair_last_abundance, rev = true) : element_isotope_pair
end

pair_last_abundance(x) = elements_abundance()[last(x)]

"""
    element_mass_delta(old_element, new_element)

Mass of `new_element` minus mass of `old_element`.
"""
element_mass_delta(old_element, new_element) = elements_mass()[new_element] - elements_mass()[old_element]

"""
    parent_mass_delta(element)

Mass of `element` minus mass of parent element of `element`.
"""
parent_mass_delta(e) = elements_mass()[e] - elements_mass()[parent_element(e)]
n_parent_mass_delta(x) = last(x) * parent_mass_delta(first(x))

"""
    deltammi(elements::Vector{Pair{String, Int}})

Sum of mass of each element in `elements` minus mass of their parent elements.
"""
deltammi(x) = sum(n_parent_mass_delta, x; init = float(0))

"""
    update_abundance(precise::Val, prev_abundance, old_element, new_element, nold, nnew, n)

Update of `prev_abundance` after `n` `old_element` replaced by `new_element`. `nold` and `nnew` are number of elements before replacing.
"""
function update_abundance(precise::Val, prev_abundance::T, old_element, new_element, nold, nnew, delta) where T
    x = get(elements_abundance(), old_element, one(T))
    y = get(elements_abundance(), new_element, one(T))
    if (x == one(T) || y == one(T))
        prev_abundance
    elseif delta == 1 
        update_abundance1(precise, prev_abundance, x, y, nold, nnew)
    else
        update_abundancen(precise, prev_abundance, x, y, nold, nnew, delta)
    end
end

update_abundance1(::Val{true}, prev_abundance, x, y, nold, nnew) = prev_abundance * (big(nold) / (nnew + 1)) * (big(y) / x)
update_abundance1(::Val{false}, prev_abundance, x, y, nold, nnew) = prev_abundance * (nold / (nnew + 1)) * (y / x)
update_abundancen(::Val{true}, prev_abundance, x, y, nold, nnew, delta) = 
    prev_abundance * (factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)) * (big(y) / x) ^ delta
function update_abundancen(::Val{false}, prev_abundance, x, y, nold, nnew, delta) 
    if (nold > nnew + delta ? check_overflow_factorial(nold, nold - delta) : check_overflow_factorial(nnew + delta, nnew))
        prev_abundance *= factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)
    else
        prev_abundance *= safe_factorial(nold, nold - delta) / safe_factorial(nnew + delta, nnew) 
    end
    convert(float(Int), prev_abundance * (y / x) ^ delta)
end

"""
    update_proportion(precise::Val, prev_proportion, nold, nnew, n) 
    
Update of `prev_proportion` after `n` elements replacing. `nold` and `nnew` are number of elements before replacing.
"""
function update_proportion(precise::Val, prev_proportion, nold, nnew, delta) 
    if delta == 1
        update_proportion1(precise, prev_proportion, nold, nnew)
    else
        update_proportionn(precise, prev_proportion, nold, nnew, delta)
    end
end

update_proportion1(::Val{true}, prev_proportion, nold, nnew) = prev_proportion * (big(nold) / (nnew + 1)) 
update_proportion1(::Val{false}, prev_proportion, nold, nnew) = prev_proportion * (nold / (nnew + 1)) 
update_proportionn(::Val{true}, prev_proportion, nold, nnew, delta) = 
    prev_proportion * (factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)) 
function update_proportionn(::Val{false}, prev_proportion, nold, nnew, delta) 
    if (nold > nnew + delta ? check_overflow_factorial(nold, nold - delta) : check_overflow_factorial(nnew + delta, nnew))
        convert(float(Int), prev_proportion * (factorial(big(nold), nold - delta) / factorial(big(nnew + delta), nnew)))
    else
        convert(float(Int), prev_proportion * (safe_factorial(nold, nold - delta) / safe_factorial(nnew + delta, nnew)))
    end
end

"""
    update_inverse_proportion(precise::Val, prev_inverse_proportion, nold, nnew, delta)

Update of `prev_inverse_proportion` after `n` elements replacing. `nold` and `nnew` are number of elements before replacing.
"""
function update_inverse_proportion(precise::Val, prev_proportion, nold, nnew, delta)
    if delta == 1 
        update_inverse_proportion1(precise, prev_proportion, nold, nnew)
    else
        update_inverse_proportionn(precise, prev_proportion, nold, nnew, delta)
    end
end

update_inverse_proportion1(::Val{true}, prev_proportion, nold, nnew) = prev_proportion * (big(nnew + 1) / nold) 
update_inverse_proportion1(::Val{false}, prev_proportion, nold, nnew) = prev_proportion * ((nnew + 1) / nold) 
update_inverse_proportionn(::Val{true}, prev_proportion, nold, nnew, delta) = 
    prev_proportion * (factorial(big(nnew + delta), nnew) / factorial(big(nold), nold - delta)) 
function update_inverse_proportionn(::Val{false}, prev_proportion, nold, nnew, delta) 
    if (nold > nnew + delta ? check_overflow_factorial(nold, nold - delta) : check_overflow_factorial(nnew + delta, nnew))
        convert(float(Int), prev_proportion * (factorial(big(nnew + delta), nnew) / factorial(big(nold), nold - delta)))
    else
        convert(float(Int), prev_proportion * (safe_factorial(nnew + delta, nnew) / safe_factorial(nold, nold - delta)))
    end
end