# ==========================================================================================================================
# Mid level MS2
function isotopologues_elements_ms2(precise::Val, it1, element_precursor_dictionary, element_product, abundance, abtype, proportion, threshold, iter, gain, loss; normalize = true)
    if isempty(element_product)
        throw(ArgumentError("Product chemical must contain at least one element."))
    elseif gain 
        msfix = get_fixmass(element_product)
        element_dictionary = get_element_dictinonary(element_precursor_dictionary)
        proportion_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * isotopicabundance(element_precursor_dictionary) / abundance
        it2 = isotopologues_elements(precise, element_dictionary, msfix, 1, Total(), proportion_cutoff)
        data = map(eachindex(it1.Element)) do i
            (; Element = it2.Element, 
            Mass = it2.Mass .+ mmi(it1.Element[i]), 
            Abundance = it2.Abundance .* (it1.Abundance[i] * proportion))
        end
    elseif hasproperty(it1, :Preab) 
        element_product_dictionary = get_element_dictinonary(element_product)
        msfix = mmi(element_product)
        element_isotope_pair = element_isotope_pairs(element_product_dictionary; sort = false)
        swap = String[]
        for (e, n) in element_product_dictionary
            m = element_precursor_dictionary[e] - n
            if m < n
                push!(swap, e)
                element_product_dictionary[e] = m
            end
        end
        @inbounds ab_dict = map(eachindex(it1.Element)) do idp
            maximal_combination_elements(precise, it1.Element[idp], element_product_dictionary, it1.Abundance[idp] * proportion * it1.Preab[idp])
        end
        proportion_cutoff = maximum(first, ab_dict) * minimum(makecrit_value(crit(threshold), abundance)) / abundance
        if iter 
            fn_main = loss ? isotopologues_elements_loss_ms2_iter : isotopologues_elements_single_ms2_iter
            fn_post = isotopologues_elements_single_ms2_post_id_iter
        else
            fn_main = isotopologues_elements_single_ms2
            fn_post = isotopologues_elements_single_ms2_post_id
        end
        fn_post_element = loss ? isotopologues_elements_single_ms2_post_loss : isotopologues_elements_single_ms2_post_product
        data = Vector{Any}(undef, length(it1.Element))
        @inbounds for (idp, el) in enumerate(it1.Element)
            element_vec, mass_vec, abundance_vec, preab_vec = fn_main(el, ab_dict[idp]..., swap, element_isotope_pair, msfix, proportion_cutoff, precise)
            element_vec, id = fn_post_element(it1, idp, element_vec, mass_vec)
            data[idp] = fn_post(element_vec, mass_vec, abundance_vec, preab_vec, id)
        end
        #         element_product_dictionary = get_element_dictinonary(element_product)
        # _, j = findmax(it1.Abundance)
        # p = maximal_combination(precise, it1.Element[j], element_product_dictionary, it1.Abundance[j] * it1.Preab[j] * proportion)
        # i = findfirst(x -> all(y -> last(y) == 0 || iselement(first(y)), x), it1.Element)
        # if !isnothing(i) 
        #     p = max(p, it1.Abundance[i] * it1.Preab[i] * proportion)
        # end
        # msfix = mmi(element_product)
        # element_isotope_pair = element_isotope_pairs(element_product_dictionary; sort = false)
        # swap = String[]
        # for (e, n) in element_product_dictionary
        #     m = element_precursor_dictionary[e] - n
        #     if m < n
        #         push!(swap, e)
        #         element_product_dictionary[e] = m
        #     end
        # end
        # # @inbounds ab_dict = map(eachindex(it1.Element)) do idp
        # #     maximal_combination_elements(precise, it1.Element[idp], element_product_dictionary, it1.Abundance[idp] * proportion * it1.Preab[idp])
        # # end
        # rcutoff = minimum(makecrit_value(crit(threshold), abundance)) / abundance
        # proportion_cutoff = p * rcutoff
        # if iter 
        #     fn_main = loss ? isotopologues_elements_loss_ms2_iter : isotopologues_elements_single_ms2_iter
        #     fn_post = isotopologues_elements_single_ms2_post_id_iter
        # else
        #     fn_main = isotopologues_elements_single_ms2
        #     fn_post = isotopologues_elements_single_ms2_post_id
        # end
        # fn_post_element = loss ? isotopologues_elements_single_ms2_post_loss : isotopologues_elements_single_ms2_post_product
        # element_product = copy(element_product_dictionary)
        # for k in keys(element_product)
        #     for e in @view elements_isotopes()[k][begin + 1:end]
        #         element_product_dictionary[e] = 0 
        #     end
        # end
        # element_product = copy(element_product_dictionary)
        # ids = sortperm(it1.Abundance; rev = true)
        # data = Vector{Any}(undef, length(ids))
        # @inbounds for idp in ids
        #     element_vec, mass_vec, abundance_vec, preab_vec = fn_main(it1.Element[idp], maximal_combination_elements!(precise, it1.Element[idp], element_product_dictionary, element_product, it1.Abundance[idp] * proportion * it1.Preab[idp])..., swap, element_isotope_pair, msfix, proportion_cutoff, precise)
        #     element_vec, id = fn_post_element(it1, idp, element_vec, mass_vec)
        #     p2 = maximum(abundance_vec)
        #     if p < p2
        #         proportion_cutoff = p2 * rcutoff
        #         p = p2
        #     end
        #     data[idp] = fn_post(element_vec, mass_vec, abundance_vec, preab_vec, id)
        # end
    else
        element_product_dictionary = get_element_dictinonary(element_product)
        msfix = mmi(element_product)
        element_isotope_pair = element_isotope_pairs(element_product_dictionary; sort = false)
        swap = String[]
        for (e, n) in element_product_dictionary
            m = element_precursor_dictionary[e] - n
            if m < n
                push!(swap, e)
                element_product_dictionary[e] = m
            end
        end
        @inbounds ab_dict = map(eachindex(it1.Element)) do idp 
            maximal_proportion_elements(precise, it1.Element[idp], element_product_dictionary, it1.Abundance[idp] * proportion)
        end
        proportion_cutoff = maximum(first, ab_dict) * minimum(makecrit_value(crit(threshold), abundance)) / abundance 
        if iter 
            fn_main = isotopologues_elements_single_ms2_iter
            fn_post = isotopologues_elements_single_ms2_post_id_iter
        else
            fn_main = isotopologues_elements_single_ms2
            fn_post = isotopologues_elements_single_ms2_post_id
        end
        data = Vector{Any}(undef, length(it1.Element))
        fn_post_element = loss ? isotopologues_elements_single_ms2_post_loss : isotopologues_elements_single_ms2_post_product
        @inbounds for (idp, el) in enumerate(it1.Element)
            element_vec, mass_vec, abundance_vec, preab_vec = fn_main(el, ab_dict[idp]..., swap, element_isotope_pair, msfix, proportion_cutoff, precise)
            element_vec, id = fn_post_element(it1, idp, element_vec, mass_vec)
            data[idp] = fn_post(element_vec, mass_vec, abundance_vec, preab_vec, id)
        end
    end
    v = getproperty.(data, :Element)
    element_product = ChainedVector(v)
    mass_product = ChainedVector(getproperty.(data, :Mass))
    abundance_pair = ChainedVector(getproperty.(data, :Abundance))
    id_pair = ChainedVector([[i for _ in eachindex(x)] for (i, x) in enumerate(v)])
    abtype = abtyped(abtype)
    # abundance_pair_normalize = normalize_abundance(abundance_pair, abundance, abtype, [Total()])
    # spectrum specific threshold ?
    if isempty(abundance_pair) 
        if iter
            (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = abundance_pair, Preab = ChainedVector(getproperty.(data, :Preab)))
        else
            (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = abundance_pair)
        end
    else
        abundance_cutoff = isempty(abundance_pair) ? 0 : minimum(makecrit_value(crit(threshold), maximum(abundance_pair)))
        id = abundance_pair .>= abundance_cutoff
        if all(id) 
            if iter 
                (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = abundance_pair, Preab = ChainedVector(getproperty.(data, :Preab)))
            else
                (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = abundance_pair) 
            end
        elseif iter
            (; ID = id_pair[id], Element = element_product[id], Mass = mass_product[id], Abundance = abundance_pair[id], Preab = ChainedVector(getproperty.(data, :Preab))[id])
        else
            (; ID = id_pair[id], Element = element_product[id], Mass = mass_product[id], Abundance = abundance_pair[id]) 
        end
    end
end

function isotopologues_elements_single_ms2(precursor_dictionary::Dict, first_abundance, element_product_dictionary::Dict, swap, element_isotope_pair, msfix, threshold, precise)
    element_vec = [swap_elements(element_product_dictionary, precursor_dictionary, swap)]
    mass_vec = [deltammi(first(element_vec)) + msfix]
    abundance_vec = [first_abundance]
    rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), threshold, swap, (true, true), precise) 
        # rec_moveisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), threshold, swap, precise)
    element_vec, mass_vec, abundance_vec, nothing
end

function isotopologues_elements_single_ms2_iter(precursor_dictionary::Dict, first_abundance, element_product_dictionary::Dict, swap, element_isotope_pair, msfix, threshold, precise)
    element_vec = [swap_elements(element_product_dictionary, precursor_dictionary, swap)]
    mass_vec = [deltammi(first(element_vec)) + msfix]
    abundance_vec = [first_abundance]
    preab_vec = [isotopologue_inverse_combination(precise, precursor_dictionary, element_product_dictionary, swap)]
    rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, swap, (true, true), precise) 
        # rec_moveisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, swap, precise)
    element_vec, mass_vec, abundance_vec, preab_vec
end

function isotopologues_elements_loss_ms2_iter(precursor_dictionary::Dict, first_abundance, element_product_dictionary::Dict, swap, element_isotope_pair, msfix, threshold, precise)
    element_vec = [swap_elements(element_product_dictionary, precursor_dictionary, swap)]
    mass_vec = [deltammi(first(element_vec)) + msfix]
    abundance_vec = [first_abundance]
    preab_vec = [isotopologue_inverse_combination(precise, precursor_dictionary, element_product_dictionary, swap)]
    rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, swap, (true, true), precise) 
        # rec_moveisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, swap, precise)
    element_vec, mass_vec, abundance_vec, preab_vec
end

function isotopologues_elements_single_ms2_post_product(it1, idp, element_vec, mass_vec) 
    id = sortperm(mass_vec)
    element_vec[id], id 
end

function isotopologues_elements_single_ms2_post_loss(it1, idp, element_vec, mass_vec) 
    id = sortperm(mass_vec)
    isotopes_precursor = it1.Isotope[idp]
    element_vec = [loss_elements(element_vec[i], isotopes_precursor) for i in id] 
    element_vec, id 
end

function isotopologues_elements_single_ms2_post_id(element_vec, mass_vec, abundance_vec, preab_vec, id)
    (; Element = element_vec, 
        Mass = mass_vec[id], 
        Abundance = abundance_vec[id]
    )
end

function isotopologues_elements_single_ms2_post_id_iter(element_vec, mass_vec, abundance_vec, preab_vec, id)
    (; Element = element_vec, 
        Mass = mass_vec[id], 
        Abundance = abundance_vec[id],
        Preab = preab_vec[id]
    )
end

# ==========================================================================================================================
# Low level MS2
function rec_exchangeisotopes!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        element_product_dictionary::Dict, 
        element_precursor_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_proportion,
        threshold,
        swap,
        current, 
        precise
    )
    next_proportion = prev_proportion
    max_proportion = prev_proportion
    (e, i) = element_isotope_pair[isotope_position]
    iter = isotope_position < lastindex(element_isotope_pair) 
    iter_e = iter && (first(element_isotope_pair[isotope_position + 1]) == e)
    if iter && (prev_proportion >= threshold || all(current))
        next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, threshold, swap, current, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    npis = get(element_product_dictionary, i, 0)
    npel = get(element_product_dictionary, e, 0)
    nrel = get(element_precursor_dictionary, e, 0) - npel
    nris = get(element_precursor_dictionary, i, 0) - npis
    if first(current)
        pis = npis 
        pel = npel 
        rel = nrel
        ris = nris
        proportion = prev_proportion
        mass = prev_mass
        backward_max_proportion = max(prev_proportion, next_proportion)
        while pis > 0 && rel >= 0
            new_proportion = update_proportion1(precise, proportion, rel, ris)
            new_proportion = update_proportion1(precise, new_proportion, pis, pel)
            if new_proportion >= threshold
                next = (true, true)
            elseif iter_e && new_proportion >= proportion
                next = max_proportion >= threshold ? (false, true) : (false, false)
            elseif new_proportion >= proportion
                next = (false, false)
            elseif iter
                max_proportion >= threshold || break
                next = (false, true)
            else
                break
            end
            pis -= 1
            pel += 1
            ris += 1 
            rel -= 1
            proportion = new_proportion
            mass += e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e)
            if !any(current)
                backward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(abundance_vec, proportion)
                push!(mass_vec, mass) 
                backward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, threshold, swap, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            elseif proportion >= threshold
                element_product_dictionary[i] = pis
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                backward_max_proportion = proportion
            elseif iter 
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                backward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, threshold, swap, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            else
                backward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, backward_max_proportion)
        end
        element_product_dictionary[i] = npis
        element_product_dictionary[e] = npel
    end
    if last(current)
        pis = npis 
        pel = npel 
        rel = nrel
        ris = nris
        proportion = prev_proportion
        mass = prev_mass
        forward_max_proportion = max(prev_proportion, next_proportion)
        while pel > 0 && ris >= 0
            new_proportion = update_proportion1(precise, proportion, pel, pis)
            new_proportion = update_proportion1(precise, new_proportion, ris, rel)
            if new_proportion >= threshold
                next = (true, true)
            elseif iter_e && new_proportion >= proportion
                next = forward_max_proportion >= threshold ? (true, false) : (false, false)
            elseif new_proportion >= proportion
                next = (false, false)
            elseif iter_e
                forward_max_proportion >= threshold || break
                next = (true, false)
            else
                break
            end
            pis += 1
            pel -= 1
            ris -= 1 
            rel += 1
            proportion = new_proportion
            mass += e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i)
            if !any(current)
                forward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, threshold, swap, next, precise)
                forward_max_proportion = max(proportion, forward_next_proportion) 
            elseif proportion >= threshold 
                element_product_dictionary[i] = pis
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_max_proportion = proportion
            elseif iter 
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                forward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, threshold, swap, next, precise) 
                forward_max_proportion = max(proportion, forward_next_proportion) 
            else
                forward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, forward_max_proportion)
        end
        element_product_dictionary[i] = npis
        element_product_dictionary[e] = npel
    end
    max_proportion
end
# function rec_exchangeisotopes!(
#         element_vec::Vector, 
#         mass_vec::Vector,
#         abundance_vec::Vector, 
#         element_product_dictionary::Dict, 
#         element_precursor_dictionary::Dict, 
#         element_isotope_pair::Vector, 
#         isotope_position::Int, 
#         prev_mass,
#         prev_proportion,
#         threshold,
#         swap,
#         current, 
#         precise
#     )
#     next_state = true
#     max_proportion = prev_proportion
#     backward_state = false
#     rec_state = false
#     backward_proportion = prev_proportion
#     (e, i) = element_isotope_pair[isotope_position]
#     if isotope_position < lastindex(element_isotope_pair) 
#         next_state, next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, threshold, swap, (true, true), precise)
#         max_proportion = max(max_proportion, next_proportion)
#     end
#     pis = get(element_product_dictionary, i, 0)
#     pel = get(element_product_dictionary, e, 0)
#     rel = get(element_precursor_dictionary, e, 0) - pel
#     ris = get(element_precursor_dictionary, i, 0) - pis
#     if first(current) && pis > 0 && rel >= 0
#         proportion = update_proportion1(precise, prev_proportion, rel, ris, 1)
#         proportion = update_proportion1(precise, proportion, pis, pel, 1)
#         backward_proportion = max(max_proportion, proportion)
#         if (proportion >= threshold) || (proportion >= prev_proportion) 
#             rec_state = true
#         elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
#             ne, ni = element_isotope_pair[isotope_position + 1]
#             # same parent replacement
#             if ne == e && get(element_precursor_dictionary, ni, 0) - get(element_product_dictionary, ni, 0) > 0
#                 rec_state = max_proportion >= threshold 
#             end
#         end
#         if rec_state 
#             new_mass = prev_mass + (e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e))
#             element_product_dictionary[i] -= 1
#             element_product_dictionary[e] += 1
#             if proportion >= threshold 
#                 push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
#                 push!(abundance_vec, proportion)
#                 push!(mass_vec, new_mass)
#             end
#             rec_state, rec_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, threshold, swap, (true, false), precise) 
#             backward_state = rec_state || backward_state 
#             backward_proportion = max(backward_proportion, rec_proportion)
#             element_product_dictionary[i] += 1
#             element_product_dictionary[e] -= 1
#         end
#     end
#     forward_state = false
#     rec_state = false
#     forward_proportion = prev_proportion
#     if last(current) && pel > 0 && ris >= 0
#         proportion = update_proportion1(precise, prev_proportion, pel, pis, 1)
#         proportion = update_proportion1(precise, proportion, ris, rel, 1)
#         forward_proportion = max(max_proportion, proportion)
#         if (proportion >= threshold) || (proportion >= prev_proportion) 
#             rec_state = true
#         elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
#             ne, ni = element_isotope_pair[isotope_position + 1]
#             # same parent replacement
#             if ne == e && get(element_product_dictionary, ni, 0) > 0
#                 rec_state = max_proportion >= threshold 
#             end
#         end
#         if rec_state
#             new_mass = prev_mass + (e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i))
#             element_product_dictionary[e] -= 1
#             element_product_dictionary[i] += 1
#             if proportion >= threshold  
#                 push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
#                 push!(abundance_vec, proportion)
#                 push!(mass_vec, new_mass)
#             end
#             rec_state, rec_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, threshold, swap, (false, true), precise) 
#             forward_state = rec_state || forward_state 
#             forward_proportion = max(forward_proportion, rec_proportion)
#             element_product_dictionary[e] += 1
#             element_product_dictionary[i] -= 1
#         end
#     end
#     forward_state || backward_state || next_state, max(forward_proportion, backward_proportion)
# end

function rec_exchangeisotopes_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_product_dictionary::Dict, 
        element_precursor_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_proportion,
        prev_preab,
        threshold,
        swap,
        current, 
        precise
    )
    next_proportion = prev_proportion
    max_proportion = prev_proportion
    (e, i) = element_isotope_pair[isotope_position]
    iter = isotope_position < lastindex(element_isotope_pair) 
    iter_e = iter && (first(element_isotope_pair[isotope_position + 1]) == e)
    if iter && (prev_proportion >= threshold || all(current))
        next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, prev_preab, threshold, swap, current, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    npis = get(element_product_dictionary, i, 0)
    npel = get(element_product_dictionary, e, 0)
    nrel = get(element_precursor_dictionary, e, 0) - npel
    nris = get(element_precursor_dictionary, i, 0) - npis
    if first(current)
        pis = npis 
        pel = npel 
        rel = nrel
        ris = nris
        preab = prev_preab
        proportion = prev_proportion
        mass = prev_mass
        backward_max_proportion = max(prev_proportion, next_proportion)
        while pis > 0 && rel >= 0
            new_proportion = update_proportion1(precise, proportion, rel, ris)
            new_proportion = update_proportion1(precise, new_proportion, pis, pel)
            if new_proportion >= threshold
                next = (true, true)
            elseif iter_e && new_proportion >= proportion
                next = max_proportion >= threshold ? (false, true) : (false, false)
            elseif new_proportion >= proportion
                next = (false, false)
            elseif iter
                max_proportion >= threshold || break
                next = (false, true)
            else
                break
            end
            preab = e in swap ? update_inverse_proportion1(precise, preab, rel, ris) : update_inverse_proportion1(precise, preab, pis, pel)
            pis -= 1
            pel += 1
            ris += 1 
            rel -= 1
            proportion = new_proportion
            mass += e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e)
            if !any(current)
                backward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass) 
                backward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, preab, threshold, swap, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            elseif proportion >= threshold
                element_product_dictionary[i] = pis
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                backward_max_proportion = proportion
            elseif iter 
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                backward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, preab, threshold, swap, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            else
                backward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, backward_max_proportion)
        end
        element_product_dictionary[i] = npis
        element_product_dictionary[e] = npel
    end
    if last(current)
        pis = npis 
        pel = npel 
        rel = nrel
        ris = nris
        preab = prev_preab
        proportion = prev_proportion
        mass = prev_mass
        forward_max_proportion = max(prev_proportion, next_proportion)
        while pel > 0 && ris >= 0
            new_proportion = update_proportion1(precise, proportion, pel, pis)
            new_proportion = update_proportion1(precise, new_proportion, ris, rel)
            if new_proportion >= threshold
                next = (true, true)
            elseif iter_e && new_proportion >= proportion
                next = forward_max_proportion >= threshold ? (true, false) : (false, false)
            elseif new_proportion >= proportion
                next = (false, false)
            elseif iter_e
                forward_max_proportion >= threshold || break
                next = (true, false)
            else
                break
            end
            preab = e in swap ? update_inverse_proportion1(precise, preab, ris, rel) : update_inverse_proportion1(precise, preab, pel, pis)
            pis += 1
            pel -= 1
            ris -= 1 
            rel += 1
            proportion = new_proportion
            mass += e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i)
            if !any(current)
                forward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, preab, threshold, swap, next, precise)
                forward_max_proportion = max(proportion, forward_next_proportion) 
            elseif proportion >= threshold 
                element_product_dictionary[i] = pis
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_max_proportion = proportion
            elseif iter 
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                forward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, threshold, preab, swap, next, precise) 
                forward_max_proportion = max(proportion, forward_next_proportion) 
            else
                forward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, forward_max_proportion)
        end
        element_product_dictionary[i] = npis
        element_product_dictionary[e] = npel
    end
    max_proportion
end

# function rec_exchangeisotopes_iter!(
#         element_vec::Vector, 
#         mass_vec::Vector,
#         abundance_vec::Vector, 
#         preab_vec::Vector, 
#         element_product_dictionary::Dict, 
#         element_precursor_dictionary::Dict, 
#         element_isotope_pair::Vector, 
#         isotope_position::Int, 
#         prev_mass,
#         prev_proportion,
#         prev_preab,
#         threshold,
#         swap,
#         backward, 
#         forward, 
#         precise
#     )
#     next_state = true
#     max_proportion = prev_proportion
#     backward_state = false
#     rec_state = false
#     backward_proportion = prev_proportion
#     (e, i) = element_isotope_pair[isotope_position]
#     if isotope_position < lastindex(element_isotope_pair) 
#         next_state, next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, prev_preab, threshold, swap, true, true, precise)
#         max_proportion = max(max_proportion, next_proportion)
#     end
#     pis = get(element_product_dictionary, i, 0)
#     pel = get(element_product_dictionary, e, 0)
#     rel = get(element_precursor_dictionary, e, 0) - pel
#     ris = get(element_precursor_dictionary, i, 0) - pis
#     if backward && pis > 0 && rel >= 0
#         proportion = update_proportion1(precise, prev_proportion, rel, ris, 1)
#         proportion = update_proportion1(precise, proportion, pis, pel, 1)
#         backward_proportion = max(max_proportion, proportion)
#         if (proportion >= threshold) || (proportion >= prev_proportion) 
#             rec_state = true
#         elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
#             ne, ni = element_isotope_pair[isotope_position + 1]
#             # same parent replacement
#             if ne == e && get(element_precursor_dictionary, ni, 0) - get(element_product_dictionary, ni, 0) > 0
#                 rec_state = max_proportion >= threshold 
#             end
#         end
#         if rec_state 
#             new_mass = prev_mass + (e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e))
#             preab = e in swap ? update_inverse_proportion1(precise, prev_preab, rel, ris, 1) : update_inverse_proportion1(precise, prev_preab, pis, pel, 1)
#             element_product_dictionary[i] -= 1
#             element_product_dictionary[e] += 1
#             if proportion >= threshold 
#                 push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
#                 push!(abundance_vec, proportion)
#                 push!(preab_vec, preab)
#                 push!(mass_vec, new_mass)
#             end
#             rec_state, rec_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, true, false, precise) 
#             backward_state = rec_state || backward_state 
#             backward_proportion = max(backward_proportion, rec_proportion)
#             element_product_dictionary[i] += 1
#             element_product_dictionary[e] -= 1
#         end
#     end
#     forward_state = false
#     rec_state = false
#     forward_proportion = prev_proportion
#     if forward && pel > 0 && ris >= 0
#         proportion = update_proportion1(precise, prev_proportion, pel, pis, 1)
#         proportion = update_proportion1(precise, proportion, ris, rel, 1)
#         forward_proportion = max(max_proportion, proportion)
#         if (proportion >= threshold) || (proportion >= prev_proportion) 
#             rec_state = true
#         elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
#             ne, ni = element_isotope_pair[isotope_position + 1]
#             # same parent replacement
#             if ne == e && get(element_product_dictionary, ni, 0) > 0
#                 rec_state = max_proportion >= threshold 
#             end
#         end
#         if rec_state
#             new_mass = prev_mass + (e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i))
#             preab = e in swap ? update_inverse_proportion1(precise, prev_preab, ris, rel, 1) : update_inverse_proportion1(precise, prev_preab, pel, pis, 1)
#             element_product_dictionary[e] -= 1
#             element_product_dictionary[i] += 1
#             if proportion >= threshold  
#                 push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
#                 push!(abundance_vec, proportion)
#                 push!(preab_vec, preab)
#                 push!(mass_vec, new_mass)
#             end
#             rec_state, rec_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, false, true, precise) 
#             forward_state = rec_state || forward_state 
#             forward_proportion = max(forward_proportion, rec_proportion)
#             element_product_dictionary[e] += 1
#             element_product_dictionary[i] -= 1
#         end
#     end
#     forward_state || backward_state || next_state, max(forward_proportion, backward_proportion)
# end

function rec_exchangeisotopes_loss_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_product_dictionary::Dict, 
        element_precursor_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_proportion,
        prev_preab,
        threshold,
        swap,
        current, 
        precise
    )
    next_proportion = prev_proportion
    max_proportion = prev_proportion
    (e, i) = element_isotope_pair[isotope_position]
    iter = isotope_position < lastindex(element_isotope_pair) 
    iter_e = iter && (first(element_isotope_pair[isotope_position + 1]) == e)
    if iter && (prev_proportion >= threshold || all(current))
        next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, prev_preab, threshold, swap, current, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    npis = get(element_product_dictionary, i, 0)
    npel = get(element_product_dictionary, e, 0)
    nrel = get(element_precursor_dictionary, e, 0) - npel
    nris = get(element_precursor_dictionary, i, 0) - npis
    if first(current)
        pis = npis 
        pel = npel 
        rel = nrel
        ris = nris
        preab = prev_preab
        proportion = prev_proportion
        mass = prev_mass
        backward_max_proportion = max(prev_proportion, next_proportion)
        while pis > 0 && rel >= 0
            new_proportion = update_proportion1(precise, proportion, rel, ris)
            new_proportion = update_proportion1(precise, new_proportion, pis, pel)
            if new_proportion >= threshold
                next = (true, true)
            elseif iter_e && new_proportion >= proportion
                next = max_proportion >= threshold ? (false, true) : (false, false)
            elseif new_proportion >= proportion
                next = (false, false)
            elseif iter
                max_proportion >= threshold || break
                next = (false, true)
            else
                break
            end
            preab = e in swap ? update_inverse_proportion1(precise, preab, rel, ris) : update_inverse_proportion1(precise, preab, pis, pel)
            pis -= 1
            pel += 1
            ris += 1 
            rel -= 1
            proportion = new_proportion
            mass += e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e)
            if !any(current)
                backward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass) 
                backward_next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, preab, threshold, swap, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            elseif proportion >= threshold
                element_product_dictionary[i] = pis
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                backward_max_proportion = proportion
            elseif iter 
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                backward_next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, preab, threshold, swap, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            else
                backward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, backward_max_proportion)
        end
        element_product_dictionary[i] = npis
        element_product_dictionary[e] = npel
    end
    if last(current)
        pis = npis 
        pel = npel 
        rel = nrel
        ris = nris
        preab = prev_preab
        proportion = prev_proportion
        mass = prev_mass
        forward_max_proportion = max(prev_proportion, next_proportion)
        while pel > 0 && ris >= 0
            new_proportion = update_proportion1(precise, proportion, pel, pis)
            new_proportion = update_proportion1(precise, new_proportion, ris, rel)
            if new_proportion >= threshold
                next = (true, true)
            elseif iter_e && new_proportion >= proportion
                next = forward_max_proportion >= threshold ? (true, false) : (false, false)
            elseif new_proportion >= proportion
                next = (false, false)
            elseif iter_e
                forward_max_proportion >= threshold || break
                next = (true, false)
            else
                break
            end
            preab = e in swap ? update_inverse_proportion1(precise, preab, ris, rel) : update_inverse_proportion1(precise, preab, pel, pis) 
            pis += 1
            pel -= 1
            ris -= 1 
            rel += 1
            proportion = new_proportion
            mass += e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i)
            if !any(current)
                forward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, preab, threshold, swap, next, precise)
                forward_max_proportion = max(proportion, forward_next_proportion) 
            elseif proportion >= threshold 
                element_product_dictionary[i] = pis
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_max_proportion = proportion
            elseif iter 
                element_product_dictionary[i] = pis
                element_product_dictionary[e] = pel
                forward_next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, mass, proportion, threshold, preab, swap, next, precise) 
                forward_max_proportion = max(proportion, forward_next_proportion) 
            else
                forward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, forward_max_proportion)
        end
        element_product_dictionary[i] = npis
        element_product_dictionary[e] = npel
    end
    max_proportion
end

# function rec_exchangeisotopes_loss_iter!(
#         element_vec::Vector, 
#         mass_vec::Vector,
#         abundance_vec::Vector, 
#         preab_vec::Vector, 
#         element_product_dictionary::Dict, 
#         element_precursor_dictionary::Dict, 
#         element_isotope_pair::Vector, 
#         isotope_position::Int, 
#         prev_mass,
#         prev_proportion,
#         prev_preab,
#         threshold,
#         swap,
#         backward, 
#         forward, 
#         precise
#     )
#     next_state = true
#     max_proportion = prev_proportion
#     backward_state = false
#     rec_state = false
#     backward_proportion = prev_proportion
#     (e, i) = element_isotope_pair[isotope_position]
#     if isotope_position < lastindex(element_isotope_pair) 
#         next_state, next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, prev_preab, threshold, swap, true, true, precise)
#         max_proportion = max(max_proportion, next_proportion)
#     end
#     pis = get(element_product_dictionary, i, 0)
#     pel = get(element_product_dictionary, e, 0)
#     rel = get(element_precursor_dictionary, e, 0) - pel
#     ris = get(element_precursor_dictionary, i, 0) - pis
#     if backward && pis > 0 && rel >= 0
#         proportion = update_proportion1(precise, prev_proportion, rel, ris, 1)
#         proportion = update_proportion1(precise, proportion, pis, pel, 1)
#         backward_proportion = max(max_proportion, proportion)
#         if (proportion >= threshold) || (proportion >= prev_proportion) 
#             rec_state = true
#         elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
#             ne, ni = element_isotope_pair[isotope_position + 1]
#             # same parent replacement
#             if ne == e 
#                 rec_state = max_proportion >= threshold 
#             end
#         end
#         if rec_state 
#             new_mass = prev_mass + (e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e))
#             preab = e in swap ? update_inverse_proportion1(precise, prev_preab, pis, pel, 1) : update_inverse_proportion1(precise, prev_preab, rel, ris, 1)
#             element_product_dictionary[i] -= 1
#             element_product_dictionary[e] += 1
#             if proportion >= threshold 
#                 push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
#                 push!(abundance_vec, proportion)
#                 push!(preab_vec, preab)
#                 push!(mass_vec, new_mass)
#             end
#             rec_state, rec_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, true, false, precise) 
#             backward_state = rec_state || backward_state 
#             backward_proportion = max(backward_proportion, rec_proportion)
#             element_product_dictionary[i] += 1
#             element_product_dictionary[e] -= 1
#         end
#     end
#     forward_state = false
#     rec_state = false
#     forward_proportion = prev_proportion
#     if forward && pel > 0 && ris >= 0
#         proportion = update_proportion1(precise, prev_proportion, pel, pis, 1)
#         proportion = update_proportion1(precise, proportion, ris, rel, 1)
#         forward_proportion = max(max_proportion, proportion)
#         if (proportion >= threshold) || (proportion >= prev_proportion) 
#             rec_state = true
#         elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
#             ne, ni = element_isotope_pair[isotope_position + 1]
#             # same parent replacement
#             if ne == e 
#                 rec_state = max_proportion >= threshold 
#             # else
#             #     # Other parent peaks
#             #     rec_state = max_proportion / prev_proportion * proportion >= threshold 
#             end
#         end
#         if rec_state
#             new_mass = prev_mass + (e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i))
#             preab = e in swap ? update_inverse_proportion1(precise, prev_preab, pel, pis, 1) : update_inverse_proportion1(precise, prev_preab, ris, rel, 1) 
#             element_product_dictionary[e] -= 1
#             element_product_dictionary[i] += 1
#             if proportion >= threshold  
#                 push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
#                 push!(abundance_vec, proportion)
#                 push!(preab_vec, preab)
#                 push!(mass_vec, new_mass)
#             end
#             rec_state, rec_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, false, true, precise) 
#             forward_state = rec_state || forward_state 
#             forward_proportion = max(forward_proportion, rec_proportion)
#             element_product_dictionary[e] += 1
#             element_product_dictionary[i] -= 1
#         end
#     end
#     forward_state || backward_state || next_state, max(forward_proportion, backward_proportion)
# end