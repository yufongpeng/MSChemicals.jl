# ==========================================================================================================================
# Mid level MS2
function isotopologues_elements_ms2(precise::Val, it1, element_precursor_dictionary, element_product, abundance, abtype, proportion, threshold, iter, gain, precursor_loss, loss, product_loss; normalize = true)
    if isempty(element_product)
        throw(ArgumentError("Product chemical must contain at least one element."))
    elseif gain 
        _, msfix = get_element_dictinonary_fixmass(element_product)
        element_dictionary, _ = get_element_dictinonary_fixmass(element_precursor_dictionary)
        proportion_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * isotopicabundance(element_precursor_dictionary) / abundance
        it2 = isotopologues_elements(precise, element_dictionary, msfix, 1, Total(), proportion_cutoff)
        data = map(eachindex(it1.Element)) do i
            (; Element = it2.Element, 
            Mass = it2.Mass .+ mmi(it1.Element[i]), 
            Abundance = it2.Abundance .* (it1.Abundance[i] * proportion))
        end
    elseif hasproperty(it1, :Preab) && (precursor_loss == loss)
        element_product_dictionary = get_element_dictinonary(element_product)
        # fp = isotopologue_combination(precise, it1.Element[j], element_product_dictionary, it1.Abundance[j] * it1.Preab[j] * proportion)
        _, j = findmax(it1.Abundance)
        p = maximal_combination(precise, it1.Element[j], element_product_dictionary, it1.Abundance[j] * it1.Preab[j] * proportion)
        i = findfirst(x -> all(iselement, keys(x)), it1.Element)
        if !isnothing(i) 
            p = max(p, it1.Abundance[i] * it1.Preab[i] * proportion)
        end
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
        proportion_cutoff = p * minimum(makecrit_value(crit(threshold), abundance)) / abundance
        if iter 
            fn_main = (loss == product_loss) ? isotopologues_elements_single_ms2_iter : isotopologues_elements_loss_ms2_iter 
            fn_post = isotopologues_elements_single_ms2_post_id_iter
        else
            fn_main = isotopologues_elements_single_ms2
            fn_post = isotopologues_elements_single_ms2_post_id
        end
        fn_post_element = loss ? isotopologues_elements_single_ms2_post_loss : isotopologues_elements_single_ms2_post_product
        @inbounds data = map(eachindex(it1.Element)) do idp
            first_abundance, first_product_dictionary = maximal_combination_elements(precise, it1.Element[idp], element_product_dictionary, it1.Abundance[idp] * proportion * it1.Preab[idp])
            element_vec, mass_vec, abundance_vec, preab_vec = fn_main(it1.Element[idp], first_product_dictionary, first_abundance, swap, element_isotope_pair, msfix, proportion_cutoff, precise)
            element_vec, id = fn_post_element(it1, idp, element_vec, mass_vec)
            fn_post(element_vec, mass_vec, abundance_vec, preab_vec, id)
        end
    else
        element_product_dictionary = get_element_dictinonary(element_product)
        _, j = findmax(it1.Abundance)
        p = maximal_proportion(precise, it1.Element[j], element_product_dictionary, it1.Abundance[j] * proportion)
        i = findfirst(x -> all(iselement, keys(x)), it1.Element)
        if !isnothing(i) 
            p = max(p, it1.Abundance[i] * proportion)
        end
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
        proportion_cutoff = p * minimum(makecrit_value(crit(threshold), abundance)) / abundance 
        # use_max = isotopologue_proportion(precise, first(it1.Element), element_product_dictionary, first(it1.Abundance) * proportion) < proportion_cutoff
        if iter 
            fn_main = isotopologues_elements_single_ms2_iter
            fn_post = isotopologues_elements_single_ms2_post_id_iter
        else
            fn_main = isotopologues_elements_single_ms2
            fn_post = isotopologues_elements_single_ms2_post_id
        end
        fn_post_element = loss ? isotopologues_elements_single_ms2_post_loss : isotopologues_elements_single_ms2_post_product
        @inbounds data = map(eachindex(it1.Element)) do idp
            first_abundance, first_product_dictionary = maximal_proportion_elements(precise, it1.Element[idp], element_product_dictionary, it1.Abundance[idp] * proportion)
            element_vec, mass_vec, abundance_vec, preab_vec = fn_main(it1.Element[idp], first_product_dictionary, first_abundance, swap, element_isotope_pair, msfix, proportion_cutoff, precise)
            element_vec, id = fn_post_element(it1, idp, element_vec, mass_vec)
            fn_post(element_vec, mass_vec, abundance_vec, preab_vec, id)
        end
    end
    v = getproperty.(data, :Element)
    element_product = ChainedVector(v)
    mass_product = ChainedVector(getproperty.(data, :Mass))
    abundance_pair = ChainedVector(getproperty.(data, :Abundance))
    id_pair = ChainedVector([[i for _ in eachindex(x)] for (i, x) in enumerate(v)])
    abtype = abtyped(abtype)
    abundance_pair_normalize = normalize_abundance(abundance_pair, abundance, abtype, [Total()])
    # spectrum specific threshold ?
    if isempty(abundance_pair_normalize) 
        if iter
            (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = normalize ? abundance_pair_normalize : abundance_pair, Preab = ChainedVector(getproperty.(data, :Preab)))
        else
            (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = normalize ? abundance_pair_normalize : abundance_pair)
        end
    else
        abundance_cutoff = isempty(abundance_pair_normalize) ? 0 : minimum(makecrit_value(crit(threshold), maximum(abundance_pair_normalize)))
        id = abundance_pair_normalize .>= abundance_cutoff
        if all(id) 
            if iter 
                (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = normalize ? abundance_pair_normalize : abundance_pair, Preab = ChainedVector(getproperty.(data, :Preab)))
            else
                (; ID = id_pair, Element = element_product, Mass = mass_product, Abundance = normalize ? abundance_pair_normalize : abundance_pair) 
            end
        elseif iter
            (; ID = id_pair[id], Element = element_product[id], Mass = mass_product[id], Abundance = normalize ? abundance_pair_normalize[id] : abundance_pair[id], Preab = ChainedVector(getproperty.(data, :Preab))[id])
        else
            (; ID = id_pair[id], Element = element_product[id], Mass = mass_product[id], Abundance = normalize ? abundance_pair_normalize[id] : abundance_pair[id]) 
        end
    end
end

function isotopologues_elements_single_ms2(precursor_dictionary::Dict, element_product_dictionary::Dict, first_abundance, swap, element_isotope_pair, msfix, threshold, precise)
    element_vec = [swap_elements(element_product_dictionary, precursor_dictionary, swap)]
    mass_vec = [deltammi(first(element_vec)) + msfix]
    abundance_vec = [first_abundance]
    rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), threshold, swap, true, true, precise) 
        # rec_moveisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), threshold, swap, precise)
    element_vec, mass_vec, abundance_vec, nothing
end

function isotopologues_elements_single_ms2_iter(precursor_dictionary::Dict, element_product_dictionary::Dict, first_abundance, swap, element_isotope_pair, msfix, threshold, precise)
    element_vec = [swap_elements(element_product_dictionary, precursor_dictionary, swap)]
    mass_vec = [deltammi(first(element_vec)) + msfix]
    abundance_vec = [first_abundance]
    preab_vec = [isotopologue_inverse_combination(precise, precursor_dictionary, element_product_dictionary, swap)]
    rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, swap, true, true, precise) 
        # rec_moveisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, swap, precise)
    element_vec, mass_vec, abundance_vec, preab_vec
end

function isotopologues_elements_loss_ms2_iter(precursor_dictionary::Dict, element_product_dictionary::Dict, first_abundance, swap, element_isotope_pair, msfix, threshold, precise)
    element_vec = [swap_elements(element_product_dictionary, precursor_dictionary, swap)]
    mass_vec = [deltammi(first(element_vec)) + msfix]
    abundance_vec = [first_abundance]
    preab_vec = [loss_inverse_combination(precise, precursor_dictionary, element_product_dictionary, swap)]
    rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, precursor_dictionary, element_isotope_pair, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, swap, true, true, precise) 
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
        backward, 
        forward, 
        precise
    )
    next_state = true
    max_proportion = prev_proportion
    backward_state = false
    rec_state = false
    backward_proportion = prev_proportion
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair) 
        next_state, next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, threshold, swap, true, true, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    pis = get!(element_product_dictionary, i, 0)
    pel = get!(element_product_dictionary, e, 0)
    rel = get(element_precursor_dictionary, e, 0) - pel
    ris = get(element_precursor_dictionary, i, 0) - pis
    if backward && pis > 0 && rel >= 0
        proportion = update_proportion(precise, prev_proportion, rel, ris, 1)
        proportion = update_proportion(precise, proportion, pis, pel, 1)
        backward_proportion = max(max_proportion, proportion)
        if (proportion >= threshold) || (proportion >= prev_proportion) 
            rec_state = true
        elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
            ne, ni = element_isotope_pair[isotope_position + 1]
            # same parent replacement
            if ne == e && get(element_precursor_dictionary, ni, 0) - get(element_product_dictionary, ni, 0) > 0
                rec_state = max_proportion >= threshold 
            end
        end
        if rec_state 
            new_mass = prev_mass + (e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e))
            element_product_dictionary[i] -= 1
            element_product_dictionary[e] += 1
            if proportion >= threshold 
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(abundance_vec, proportion)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, threshold, swap, true, false, precise) 
            backward_state = rec_state || backward_state 
            backward_proportion = max(backward_proportion, rec_proportion)
            element_product_dictionary[i] += 1
            element_product_dictionary[e] -= 1
        end
    end
    forward_state = false
    rec_state = false
    forward_proportion = prev_proportion
    if forward && pel > 0 && ris >= 0
        proportion = update_proportion(precise, prev_proportion, pel, pis, 1)
        proportion = update_proportion(precise, proportion, ris, rel, 1)
        forward_proportion = max(max_proportion, proportion)
        if (proportion >= threshold) || (proportion >= prev_proportion) 
            rec_state = true
        elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
            ne, ni = element_isotope_pair[isotope_position + 1]
            # same parent replacement
            if ne == e && get(element_product_dictionary, ni, 0) > 0
                rec_state = max_proportion >= threshold 
            end
        end
        if rec_state
            new_mass = prev_mass + (e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i))
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
            if proportion >= threshold  
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(abundance_vec, proportion)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, threshold, swap, false, true, precise) 
            forward_state = rec_state || forward_state 
            forward_proportion = max(forward_proportion, rec_proportion)
            element_product_dictionary[e] += 1
            element_product_dictionary[i] -= 1
        end
    end
    forward_state || backward_state || next_state, max(forward_proportion, backward_proportion)
end

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
        backward, 
        forward, 
        precise
    )
    next_state = true
    max_proportion = prev_proportion
    backward_state = false
    rec_state = false
    backward_proportion = prev_proportion
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair) 
        next_state, next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, prev_preab, threshold, swap, true, true, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    pis = get!(element_product_dictionary, i, 0)
    pel = get!(element_product_dictionary, e, 0)
    rel = get(element_precursor_dictionary, e, 0) - pel
    ris = get(element_precursor_dictionary, i, 0) - pis
    if backward && pis > 0 && rel >= 0
        proportion = update_proportion(precise, prev_proportion, rel, ris, 1)
        proportion = update_proportion(precise, proportion, pis, pel, 1)
        backward_proportion = max(max_proportion, proportion)
        if (proportion >= threshold) || (proportion >= prev_proportion) 
            rec_state = true
        elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
            ne, ni = element_isotope_pair[isotope_position + 1]
            # same parent replacement
            if ne == e && get(element_precursor_dictionary, ni, 0) - get(element_product_dictionary, ni, 0) > 0
                rec_state = max_proportion >= threshold 
            end
        end
        if rec_state 
            new_mass = prev_mass + (e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e))
            preab = e in swap ? update_inverse_proportion(precise, prev_preab, rel, ris, 1) : update_inverse_proportion(precise, prev_preab, pis, pel, 1)
            element_product_dictionary[i] -= 1
            element_product_dictionary[e] += 1
            if proportion >= threshold 
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(abundance_vec, proportion)
                push!(preab_vec, preab)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, true, false, precise) 
            backward_state = rec_state || backward_state 
            backward_proportion = max(backward_proportion, rec_proportion)
            element_product_dictionary[i] += 1
            element_product_dictionary[e] -= 1
        end
    end
    forward_state = false
    rec_state = false
    forward_proportion = prev_proportion
    if forward && pel > 0 && ris >= 0
        proportion = update_proportion(precise, prev_proportion, pel, pis, 1)
        proportion = update_proportion(precise, proportion, ris, rel, 1)
        forward_proportion = max(max_proportion, proportion)
        if (proportion >= threshold) || (proportion >= prev_proportion) 
            rec_state = true
        elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
            ne, ni = element_isotope_pair[isotope_position + 1]
            # same parent replacement
            if ne == e && get(element_product_dictionary, ni, 0) > 0
                rec_state = max_proportion >= threshold 
            end
        end
        if rec_state
            new_mass = prev_mass + (e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i))
            preab = e in swap ? update_inverse_proportion(precise, prev_preab, ris, rel, 1) : update_inverse_proportion(precise, prev_preab, pel, pis, 1)
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
            if proportion >= threshold  
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(abundance_vec, proportion)
                push!(preab_vec, preab)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, false, true, precise) 
            forward_state = rec_state || forward_state 
            forward_proportion = max(forward_proportion, rec_proportion)
            element_product_dictionary[e] += 1
            element_product_dictionary[i] -= 1
        end
    end
    forward_state || backward_state || next_state, max(forward_proportion, backward_proportion)
end

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
        backward, 
        forward, 
        precise
    )
    next_state = true
    max_proportion = prev_proportion
    backward_state = false
    rec_state = false
    backward_proportion = prev_proportion
    (e, i) = element_isotope_pair[isotope_position]
    if isotope_position < lastindex(element_isotope_pair) 
        next_state, next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_proportion, prev_preab, threshold, swap, true, true, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    pis = get!(element_product_dictionary, i, 0)
    pel = get!(element_product_dictionary, e, 0)
    rel = get(element_precursor_dictionary, e, 0) - pel
    ris = get(element_precursor_dictionary, i, 0) - pis
    if backward && pis > 0 && rel >= 0
        proportion = update_proportion(precise, prev_proportion, rel, ris, 1)
        proportion = update_proportion(precise, proportion, pis, pel, 1)
        backward_proportion = max(max_proportion, proportion)
        if (proportion >= threshold) || (proportion >= prev_proportion) 
            rec_state = true
        elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
            ne, ni = element_isotope_pair[isotope_position + 1]
            # same parent replacement
            if ne == e && get(element_precursor_dictionary, ni, 0) - get(element_product_dictionary, ni, 0) > 0
                rec_state = max_proportion >= threshold 
            end
        end
        if rec_state 
            new_mass = prev_mass + (e in swap ? element_mass_delta(e, i) : element_mass_delta(i, e))
            preab = e in swap ? update_inverse_proportion(precise, prev_preab, pis, pel, 1) : update_inverse_proportion(precise, prev_preab, rel, ris, 1)
            element_product_dictionary[i] -= 1
            element_product_dictionary[e] += 1
            if proportion >= threshold 
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(abundance_vec, proportion)
                push!(preab_vec, preab)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, true, false, precise) 
            backward_state = rec_state || backward_state 
            backward_proportion = max(backward_proportion, rec_proportion)
            element_product_dictionary[i] += 1
            element_product_dictionary[e] -= 1
        end
    end
    forward_state = false
    rec_state = false
    forward_proportion = prev_proportion
    if forward && pel > 0 && ris >= 0
        proportion = update_proportion(precise, prev_proportion, pel, pis, 1)
        proportion = update_proportion(precise, proportion, ris, rel, 1)
        forward_proportion = max(max_proportion, proportion)
        if (proportion >= threshold) || (proportion >= prev_proportion) 
            rec_state = true
        elseif (prev_proportion >= proportion > 0) && isotope_position < lastindex(element_isotope_pair) 
            ne, ni = element_isotope_pair[isotope_position + 1]
            # same parent replacement
            if ne == e && get(element_product_dictionary, ni, 0) > 0
                rec_state = max_proportion >= threshold 
            # else
            #     # Other parent peaks
            #     rec_state = max_proportion / prev_proportion * proportion >= threshold 
            end
        end
        if rec_state
            new_mass = prev_mass + (e in swap ? element_mass_delta(i, e) : element_mass_delta(e, i))
            preab = e in swap ? update_inverse_proportion(precise, prev_preab, pel, pis, 1) : update_inverse_proportion(precise, prev_preab, ris, rel, 1) 
            element_product_dictionary[e] -= 1
            element_product_dictionary[i] += 1
            if proportion >= threshold  
                push!(element_vec, swap_elements(element_product_dictionary, element_precursor_dictionary, swap))
                push!(abundance_vec, proportion)
                push!(preab_vec, preab)
                push!(mass_vec, new_mass)
            end
            rec_state, rec_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product_dictionary, element_precursor_dictionary, element_isotope_pair, isotope_position, new_mass, proportion, preab, threshold, swap, false, true, precise) 
            forward_state = rec_state || forward_state 
            forward_proportion = max(forward_proportion, rec_proportion)
            element_product_dictionary[e] += 1
            element_product_dictionary[i] -= 1
        end
    end
    forward_state || backward_state || next_state, max(forward_proportion, backward_proportion)
end