# ==========================================================================================================================
# Mid level MS2
function isotopologues_elements_ms2(precise::Val, it1, element_precursor_dictionary, element_precursor, element_product, abundance, abtype, proportion, threshold, iter, gain, loss; normalize = true)
    if isempty(element_product)
        throw(ArgumentError("Product chemical must contain at least one element."))
    elseif gain 
        msfix = get_fixmass(element_product)
        element_dictionary = get_element_dictionary(element_precursor)
        proportion_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * isotopicabundance(element_precursor) / abundance
        it2 = isotopologues_elements_ms1(precise, element_dictionary, msfix, 1, Total(), proportion_cutoff)
        data = map(eachindex(it1.Element)) do i
            (; Element = it2.Element, 
            Mass = it2.Mass .+ mmi(it1.Element[i]), 
            Abundance = it2.Abundance .* (it1.Abundance[i] * proportion))
        end
    elseif hasproperty(it1, :Preab) 
        element_product_vp = get_element(element_product)
        nisotopes = get_nisotopes(element_product_vp)
        isotopes = get_isotopes(element_product_vp)
        # for (i, k) in enumerate(isotopes) 
        #     ismajor(k) || continue 
        #     isotopes[i] = parent_element(k)
        # end
        iid = findall(!ismajor, isotopes)
        msfix = mmi(element_product)
        leftover = [k for (k, v) in element_precursor_dictionary if isnothing(findfirst(x -> first(x) == k, element_product_vp))]
        # element_isotope_pair = element_isotope_pairs(element_product_dictionary; sort = false)
        @inbounds ab_dict = map(eachindex(it1.Element)) do idp
            maximal_combination_composition(precise, it1.Element[idp], element_product_vp, leftover, it1.Abundance[idp] * proportion * it1.Preab[idp])
        end
        proportion_cutoff = maximum(first, ab_dict) * minimum(makecrit_value(crit(threshold), abundance)) / abundance
        if iter 
            fn_main = loss ? isotopologues_elements_loss_ms2_iter : isotopologues_elements_single_ms2_iter
            fn_post = isotopologues_elements_single_ms2_post_id_iter
        else
            fn_main = isotopologues_elements_single_ms2
            fn_post = isotopologues_elements_single_ms2_post_id
        end
        detected_isotopes = isotopes[iid]
        if loss 
            fn_post_element = isotopologues_elements_single_ms2_post_loss 
            precursor_isotopes = filter!(!ismajor, get_isotopes(collect(element_precursor_dictionary)))
            product_isotopes = precursor_isotopes
        else
            fn_post_element = isotopologues_elements_single_ms2_post_product
            precursor_isotopes = nothing
            product_isotopes = detected_isotopes
        end
        data = Vector{Any}(undef, length(it1.Element))
        @inbounds for (idp, el) in enumerate(it1.Element)
            element_vec, mass_vec, abundance_vec, preab_vec = fn_main(el, ab_dict[idp]..., isotopes, nisotopes, iid, msfix, proportion_cutoff, precise)
            element_vec, id = fn_post_element(it1, idp, element_vec, mass_vec, detected_isotopes, precursor_isotopes)
            data[idp] = fn_post(element_vec, mass_vec, abundance_vec, preab_vec, id)
        end
    else
        element_product_vp = get_element(element_product)
        nisotopes = get_nisotopes(element_product_vp)
        isotopes = get_isotopes(element_product_vp)
        # for (i, k) in enumerate(isotopes) 
        #     ismajor(k) || continue 
        #     isotopes[i] = parent_element(k)
        # end
        iid = findall(!ismajor, isotopes)
        msfix = mmi(element_product)
        leftover = [k for (k, v) in element_precursor_dictionary if isnothing(findfirst(x -> first(x) == k, element_product_vp))]
        # element_isotope_pair = element_isotope_pairs(element_product_vp; sort = false)
        @inbounds ab_dict = map(eachindex(it1.Element)) do idp 
            maximal_proportion_composition(precise, it1.Element[idp], element_product_vp, leftover, it1.Abundance[idp] * proportion)
        end
        proportion_cutoff = maximum(first, ab_dict) * minimum(makecrit_value(crit(threshold), abundance)) / abundance 
        if iter 
            fn_main = isotopologues_elements_single_ms2_iter
            fn_post = isotopologues_elements_single_ms2_post_id_iter
        else
            fn_main = isotopologues_elements_single_ms2
            fn_post = isotopologues_elements_single_ms2_post_id
        end
        detected_isotopes = isotopes[iid]
        if loss 
            fn_post_element = isotopologues_elements_single_ms2_post_loss 
            precursor_isotopes = filter!(!ismajor, get_isotopes(collect(element_precursor_dictionary)))
            product_isotopes = precursor_isotopes
        else
            fn_post_element = isotopologues_elements_single_ms2_post_product
            precursor_isotopes = nothing
            product_isotopes = detected_isotopes
        end
        data = Vector{Any}(undef, length(it1.Element))
        @inbounds for (idp, el) in enumerate(it1.Element)
            element_vec, mass_vec, abundance_vec, preab_vec = fn_main(el, ab_dict[idp]..., isotopes, nisotopes, iid, msfix, proportion_cutoff, precise)
            element_vec, id = fn_post_element(it1, idp, element_vec, mass_vec, detected_isotopes, precursor_isotopes)
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
            (; ID = id_pair, Element = ElementsVector[], Mass = mass_product, Abundance = abundance_pair, Preab = ChainedVector(getproperty.(data, :Preab)))
        else
            (; ID = id_pair, Element = ElementsVector[], Mass = mass_product, Abundance = abundance_pair)
        end
    else
        abundance_cutoff = isempty(abundance_pair) ? 0 : minimum(makecrit_value(crit(threshold), maximum(abundance_pair)))
        id = abundance_pair .>= abundance_cutoff
        if all(id) 
            if iter 
                (; ID = id_pair, Element = gain ? element_product : [ElementsVector(product_isotopes, x) for x in element_product], Mass = mass_product, Abundance = abundance_pair, Preab = ChainedVector(getproperty.(data, :Preab)))
            else
                (; ID = id_pair, Element = gain ? element_product : [ElementsVector(product_isotopes, x) for x in element_product], Mass = mass_product, Abundance = abundance_pair) 
            end
        elseif iter
            (; ID = id_pair[id], Element = gain ? element_product[id] : [ElementsVector(product_isotopes, x) for x in element_product[id]], Mass = mass_product[id], Abundance = abundance_pair[id], Preab = ChainedVector(getproperty.(data, :Preab))[id])
        else
            (; ID = id_pair[id], Element = gain ? element_product[id] : [ElementsVector(product_isotopes, x) for x in element_product[id]], Mass = mass_product[id], Abundance = abundance_pair[id]) 
        end
    end
end

function isotopologues_elements_single_ms2(precursor_dictionary::Dict, first_abundance, element_product, isotopes, nisotopes, iid, msfix, threshold, precise)
    element_vec = [element_product[iid]]
    mass_vec = [deltammi(isotopes, first(element_vec)) + msfix]
    abundance_vec = [first_abundance]
    rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product, precursor_dictionary, isotopes, nisotopes, iid, 2, 1, 1, first(mass_vec), first(abundance_vec), threshold, (true, true), precise) 
    element_vec, mass_vec, abundance_vec, nothing
end

function isotopologues_elements_single_ms2_iter(precursor_dictionary::Dict, first_abundance, element_product, isotopes, nisotopes, iid, msfix, threshold, precise)
    element_vec = [element_product[iid]]
    mass_vec = [deltammi(isotopes, first(element_vec)) + msfix]
    abundance_vec = [first_abundance]
    preab_vec = [isotopologue_inverse_combination(precise, element_product, isotopes)]
    rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product, precursor_dictionary, isotopes, nisotopes, iid, 2, 1, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, (true, true), precise) 
    element_vec, mass_vec, abundance_vec, preab_vec
end

function isotopologues_elements_loss_ms2_iter(precursor_dictionary::Dict, first_abundance, element_product, isotopes, nisotopes, iid, msfix, threshold, precise)
    element_vec = [element_product[iid]]
    mass_vec = [deltammi(isotopes, first(element_vec)) + msfix]
    abundance_vec = [first_abundance]
    preab_vec = [isotopologue_inverse_combination(precise, element_product, isotopes)]
    rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product, precursor_dictionary, isotopes, nisotopes, iid, 2, 1, 1, first(mass_vec), first(abundance_vec), first(preab_vec), threshold, (true, true), precise) 
    element_vec, mass_vec, abundance_vec, preab_vec
end

function isotopologues_elements_single_ms2_post_product(it1, idp, element_vec, mass_vec, isotopes, precursor_isotopes) 
    id = sortperm(mass_vec)
    element_vec[id], id 
end

function isotopologues_elements_single_ms2_post_loss(it1, idp, element_vec, mass_vec, isotopes, precursor_isotopes) 
    id = sortperm(mass_vec)
    precursor_vec = zeros(Int, length(precursor_isotopes))
    for (k, v) in it1.Element[idp]
        i = findfirst(==(k), precursor_isotopes)
        isnothing(i) && continue
        precursor_vec[i] = -v 
    end
    element_vec = number_loss_elements(isotopes, element_vec, precursor_isotopes, precursor_vec, id) 
    element_vec, id 
end

function number_loss_elements(product_isotopes, product_vec, precursor_isotopes, precursor_vec, id)
    new_vps = [copy(precursor_vec) for _ in id]
    for (ie, e) in enumerate(product_isotopes)
        i = findfirst(==(e), precursor_isotopes)
        isnothing(i) && continue
        for (new_vp, j) in zip(new_vps, id)
            new_vp[i] += product_vec[j][ie]
        end
    end
    new_vps
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
        element_vp::Vector, 
        element_precursor_dictionary::Dict, 
        isotopes::Vector,
        nisotopes::Vector, 
        iid::Vector, 
        isotope_position::Int, 
        element_position::Int,
        n_position::Int,
        prev_mass,
        prev_proportion,
        threshold,
        current, 
        precise
    )
    if nisotopes[n_position] <= isotope_position - element_position
        element_position = isotope_position
        n_position += 1
        isotope_position = element_position + 1
    end
    while nisotopes[n_position] == 1 
        n_position += 1
        n_position > lastindex(nisotopes) && return prev_proportion
        element_position += 1
        isotope_position += 1
    end
    isotope_position > lastindex(isotopes) && return prev_proportion
    next_proportion = prev_proportion
    max_proportion = prev_proportion
    iter = isotope_position < lastindex(isotopes)
    iter_e = iter && (isotope_position < element_position + nisotopes[n_position] - 1)
    if iter && (prev_proportion >= threshold || all(current))
        next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, prev_mass, prev_proportion, threshold, current, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    e = isotopes[element_position]
    i = isotopes[isotope_position]
    npis = element_vp[isotope_position]
    npel = element_vp[element_position]
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
            elseif new_proportion >= proportion
                next = (iter_e && backward_max_proportion >= threshold) ? (false, true) : (false, false)
            elseif iter
                backward_max_proportion >= threshold || break
                next = (false, true)
            else
                break
            end
            pis -= 1
            pel += 1
            ris += 1 
            rel -= 1
            proportion = new_proportion
            mass += element_mass_delta(i, e)
            if !any(current)
                backward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                push!(element_vec, element_vp[iid])
                push!(abundance_vec, proportion)
                push!(mass_vec, mass) 
                backward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, threshold, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            elseif proportion >= threshold
                element_vp[isotope_position] = pis
                push!(element_vec, element_vp[iid])
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                backward_max_proportion = proportion
            elseif iter 
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                backward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, threshold, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            else
                backward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, backward_max_proportion)
        end
        element_vp[isotope_position] = npis
        element_vp[element_position] = npel
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
            elseif new_proportion >= proportion
                next = (iter_e && forward_max_proportion >= threshold) ? (true, false) : (false, false)
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
            mass += element_mass_delta(e, i)
            if !any(current)
                forward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                push!(element_vec, element_vp[iid])
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, threshold, next, precise)
                forward_max_proportion = max(proportion, forward_next_proportion) 
            elseif proportion >= threshold 
                element_vp[isotope_position] = pis
                push!(element_vec, element_vp[iid])
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_max_proportion = proportion
            elseif iter 
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                forward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, threshold, next, precise) 
                forward_max_proportion = max(proportion, forward_next_proportion) 
            else
                forward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, forward_max_proportion)
        end
        element_vp[isotope_position] = npis
        element_vp[element_position] = npel
    end
    max_proportion
end

function rec_exchangeisotopes_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_vp::Vector, 
        element_precursor_dictionary::Dict, 
        isotopes::Vector,
        nisotopes::Vector, 
        iid::Vector, 
        isotope_position::Int, 
        element_position::Int,
        n_position::Int,
        prev_mass,
        prev_proportion,
        prev_preab,
        threshold,
        current, 
        precise
    )
    if nisotopes[n_position] <= isotope_position - element_position
        element_position = isotope_position
        n_position += 1
        isotope_position = element_position + 1
    end
    while nisotopes[n_position] == 1 
        n_position += 1
        n_position > lastindex(nisotopes) && return prev_proportion
        element_position += 1
        isotope_position += 1
    end
    isotope_position > lastindex(isotopes) && return prev_proportion
    next_proportion = prev_proportion
    max_proportion = prev_proportion
    iter = isotope_position < lastindex(isotopes)
    iter_e = iter && (isotope_position < element_position + nisotopes[n_position] - 1)
    if iter && (prev_proportion >= threshold || all(current))
        next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, prev_mass, prev_proportion, prev_preab, threshold, current, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    e = isotopes[element_position]
    i = isotopes[isotope_position]
    npis = element_vp[isotope_position]
    npel = element_vp[element_position]
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
            elseif new_proportion >= proportion
                next = (iter_e && backward_max_proportion >= threshold) ? (false, true) : (false, false)
            elseif iter
                backward_max_proportion >= threshold || break
                next = (false, true)
            else
                break
            end
            preab = update_inverse_proportion1(precise, preab, pis, pel)
            pis -= 1
            pel += 1
            ris += 1 
            rel -= 1
            proportion = new_proportion
            mass += element_mass_delta(i, e)
            if !any(current)
                backward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass) 
                backward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, preab, threshold, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            elseif proportion >= threshold
                element_vp[isotope_position] = pis
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                backward_max_proportion = proportion
            elseif iter 
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                backward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, preab, threshold, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            else
                backward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, backward_max_proportion)
        end
        element_vp[isotope_position] = npis
        element_vp[element_position] = npel
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
            elseif new_proportion >= proportion
                next = (iter_e && forward_max_proportion >= threshold) ? (true, false) : (false, false)
            elseif iter_e
                forward_max_proportion >= threshold || break
                next = (true, false)
            else
                break
            end
            preab = update_inverse_proportion1(precise, preab, pel, pis)
            pis += 1
            pel -= 1
            ris -= 1 
            rel += 1
            proportion = new_proportion
            mass += element_mass_delta(e, i)
            if !any(current)
                forward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, preab, threshold, next, precise)
                forward_max_proportion = max(proportion, forward_next_proportion) 
            elseif proportion >= threshold 
                element_vp[isotope_position] = pis
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_max_proportion = proportion
            elseif iter 
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                forward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, threshold, preab, next, precise) 
                forward_max_proportion = max(proportion, forward_next_proportion) 
            else
                forward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, forward_max_proportion)
        end
        element_vp[isotope_position] = npis
        element_vp[element_position] = npel
    end
    max_proportion
end

function rec_exchangeisotopes_loss_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_vp::Vector, 
        element_precursor_dictionary::Dict, 
        isotopes::Vector,
        nisotopes::Vector, 
        iid::Vector, 
        isotope_position::Int, 
        element_position::Int,
        n_position::Int,
        prev_mass,
        prev_proportion,
        prev_preab,
        threshold,
        current, 
        precise
    )
    if nisotopes[n_position] <= isotope_position - element_position
        element_position = isotope_position
        n_position += 1
        isotope_position = element_position + 1
    end
    while nisotopes[n_position] == 1 
        n_position += 1
        n_position > lastindex(nisotopes) && return prev_proportion
        element_position += 1
        isotope_position += 1
    end
    isotope_position > lastindex(isotopes) && return prev_proportion
    next_proportion = prev_proportion
    max_proportion = prev_proportion
    iter = isotope_position < lastindex(isotopes)
    iter_e = iter && (isotope_position < element_position + nisotopes[n_position] - 1)
    if iter && (prev_proportion >= threshold || all(current))
        next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, prev_mass, prev_proportion, prev_preab, threshold, current, precise)
        max_proportion = max(max_proportion, next_proportion)
    end
    e = isotopes[element_position]
    i = isotopes[isotope_position]
    npis = element_vp[isotope_position]
    npel = element_vp[element_position]
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
            elseif new_proportion >= proportion
                next = (iter_e && backward_max_proportion >= threshold) ? (false, true) : (false, false)
            elseif iter
                backward_max_proportion >= threshold || break
                next = (false, true)
            else
                break
            end
            preab = update_inverse_proportion1(precise, preab, pis, pel)
            pis -= 1
            pel += 1
            ris += 1 
            rel -= 1
            proportion = new_proportion
            mass += element_mass_delta(i, e)
            if !any(current)
                backward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass) 
                backward_next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, preab, threshold, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            elseif proportion >= threshold
                element_vp[isotope_position] = pis
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                backward_max_proportion = proportion
            elseif iter 
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                backward_next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, preab, threshold, next, precise)
                backward_max_proportion = max(proportion, backward_next_proportion)
            else
                backward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, backward_max_proportion)
        end
        element_vp[isotope_position] = npis
        element_vp[element_position] = npel
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
            elseif new_proportion >= proportion
                next = (iter_e && forward_max_proportion >= threshold) ? (true, false) : (false, false)
            elseif iter_e
                forward_max_proportion >= threshold || break
                next = (true, false)
            else
                break
            end
            preab = update_inverse_proportion1(precise, preab, pel, pis) 
            pis += 1
            pel -= 1
            ris -= 1 
            rel += 1
            proportion = new_proportion
            mass += element_mass_delta(e, i)
            if !any(current)
                forward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, preab, threshold, next, precise)
                forward_max_proportion = max(proportion, forward_next_proportion) 
            elseif proportion >= threshold 
                element_vp[isotope_position] = pis
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_max_proportion = proportion
            elseif iter 
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                forward_next_proportion = rec_exchangeisotopes_loss_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, threshold, preab, next, precise) 
                forward_max_proportion = max(proportion, forward_next_proportion) 
            else
                forward_max_proportion = proportion
            end
            max_proportion = max(max_proportion, forward_max_proportion)
        end
        element_vp[isotope_position] = npis
        element_vp[element_position] = npel
    end
    max_proportion
end