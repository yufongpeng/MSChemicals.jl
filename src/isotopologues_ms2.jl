# ==========================================================================================================================
# Mid level MS2
function isotopologues_elements_ms2(precise::Val, precursor, product_sch, product, it1, element_precursor_dictionary, element_precursor, element_product, abundance, abtype, proportion, threshold, iter, gain, loss, sort; normalize = true)
    net_charge = charge(product)
    abs_charge = max(1, abs(net_charge))
    if isempty(element_product)
        throw(ArgumentError("Product chemical must contain at least one element."))
    elseif gain 
        msfix = get_fixmass(element_product)
        element_dictionary = get_element_dictionary(element_precursor)
        proportion_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * isotopicabundance(element_precursor) / abundance
        it2 = isotopologues_elements_ms1(precise, product_sch, element_dictionary, msfix, 1, Total(), proportion_cutoff, iter, sort; net_charge)
        abundance_cutoff = minimum(makecrit_value(crit(threshold), abundance)) * maximum(it2.Abundance) * maximum(it1.Abundance) / abundance
        data = [
            Vector{Vector{<:Vector}}(undef, length(precursor)),
            Vector{Vector{float(Int)}}(undef, length(precursor)),
            Vector{Vector{typeof(return_abundance(precise, big(0.0)))}}(undef, length(precursor))
        ]
        for idp in eachindex(IndexLinear(), precursor)
            d1 = Vector[]
            d2 = float(Int)[]
            d3 = typeof(return_abundance(precise, big(0.0)))[]
            msf = mmi(it1.Element[idp]) / abs_charge
            abf = it1.Abundance[idp] * proportion
            for (i, ab) in pairs(IndexLinear(), it2.Abundance)
                abn = ab * abf 
                abn < abundance_cutoff && continue 
                push!(d1, vcat(precursor[idp], it2.Chemical[i]))
                push!(d2, it2.Mass[i] + msf)
                push!(d3, abn)
            end
            data[1][idp] = d1
            data[2][idp] = d2
            data[3][idp] = d3
        end
    elseif hasproperty(it1, :Preab) 
        element_product_vp = get_element(element_product)
        nisotopes = get_nisotopes(element_product_vp)
        isotopes = get_isotopes(element_product_vp)
        iid = findall(!ismajor, isotopes)
        msfix = mmi(element_product) - (net_charge == 0 ? 0 : net_charge * ME)
        leftover = [k for (k, v) in element_precursor_dictionary if isnothing(findfirst(x -> first(x) == k, element_product_vp))]
        @inbounds ab_dict = map(eachindex(IndexLinear(), precursor)) do idp
            maximal_combination_composition(precise, it1.Element[idp], element_product_vp, leftover, it1.Abundance[idp] * proportion * it1.Preab[idp])
        end
        proportion_cutoff = maximum(first, ab_dict) * minimum(makecrit_value(crit(threshold), abundance)) / abundance
        if iter 
            fn_main = isotopologues_elements_single_ms2_iter
            fn_push = sort ? push_data_iter_sort! : push_data_iter!
            data = [
                Vector{Vector{<:Vector}}(undef, length(precursor)),
                Vector{Vector{float(Int)}}(undef, length(precursor)),
                Vector{Vector{typeof(return_abundance(precise, big(0.0)))}}(undef, length(precursor)),
                Vector{Vector{typeof(return_abundance(precise, big(0.0)))}}(undef, length(precursor)),
            ]
        else
            fn_main = isotopologues_elements_single_ms2
            fn_push = sort ? push_data_sort! : push_data!
            data = [
                Vector{Vector{<:Vector}}(undef, length(precursor)),
                Vector{Vector{float(Int)}}(undef, length(precursor)),
                Vector{Vector{typeof(return_abundance(precise, big(0.0)))}}(undef, length(precursor))
            ]
        end
        detected_isotopes = isotopes[iid]
        if loss 
            fn_id = sort ? loss_id_sort : loss_id
            precursor_isotopes = filter!(!ismajor, get_isotopes(collect(element_precursor_dictionary)))
            product_isotopes = precursor_isotopes
        else
            fn_id = sort ? product_id_sort : product_id
            precursor_isotopes = nothing
            product_isotopes = detected_isotopes
        end
        @inbounds for idp in eachindex(IndexLinear(), precursor)
            element_vec, mass_vec, abundance_vec, preab_vec = fn_main(it1.Element[idp], ab_dict[idp]..., isotopes, nisotopes, iid, msfix, proportion_cutoff, precise, abs_charge)
            chemical_vec, id = fn_id(it1, idp, element_vec, mass_vec, detected_isotopes, precursor_isotopes, product_isotopes, precursor, product_sch)
            fn_push(data, idp, chemical_vec, mass_vec, abundance_vec, preab_vec, id)
        end
    else
        element_product_vp = get_element(element_product)
        nisotopes = get_nisotopes(element_product_vp)
        isotopes = get_isotopes(element_product_vp)
        iid = findall(!ismajor, isotopes)
        msfix = mmi(element_product) - (net_charge == 0 ? 0 : net_charge * ME)
        leftover = [k for (k, v) in element_precursor_dictionary if isnothing(findfirst(x -> first(x) == k, element_product_vp))]
        @inbounds ab_dict = map(eachindex(IndexLinear(), precursor)) do idp 
            maximal_proportion_composition(precise, it1.Element[idp], element_product_vp, leftover, it1.Abundance[idp] * proportion)
        end
        proportion_cutoff = maximum(first, ab_dict) * minimum(makecrit_value(crit(threshold), abundance)) / abundance 
        if iter 
            fn_main = isotopologues_elements_single_ms2_iter
            fn_push = sort ? push_data_iter_sort! : push_data_iter!
            data = [
                Vector{Vector{<:Vector}}(undef, length(precursor)),
                Vector{Vector{float(Int)}}(undef, length(precursor)),
                Vector{Vector{typeof(return_abundance(precise, big(0.0)))}}(undef, length(precursor)),
                Vector{Vector{typeof(return_abundance(precise, big(0.0)))}}(undef, length(precursor)),
            ]
        else
            fn_main = isotopologues_elements_single_ms2
            fn_push = sort ? push_data_sort! : push_data!
            data = [
                Vector{Vector{<:Vector}}(undef, length(precursor)),
                Vector{Vector{float(Int)}}(undef, length(precursor)),
                Vector{Vector{typeof(return_abundance(precise, big(0.0)))}}(undef, length(precursor))
            ]
        end
        detected_isotopes = isotopes[iid]
        if loss 
            fn_id = sort ? loss_id_sort : loss_id
            precursor_isotopes = filter!(!ismajor, get_isotopes(collect(element_precursor_dictionary)))
            product_isotopes = precursor_isotopes
        else
            fn_id = sort ? product_id_sort : product_id
            precursor_isotopes = nothing
            product_isotopes = detected_isotopes
        end
        @inbounds for idp in eachindex(IndexLinear(), precursor)
            element_vec, mass_vec, abundance_vec, preab_vec = fn_main(it1.Element[idp], ab_dict[idp]..., isotopes, nisotopes, iid, msfix, proportion_cutoff, precise, abs_charge)
            chemical_vec, id = fn_id(it1, idp, element_vec, mass_vec, detected_isotopes, precursor_isotopes, product_isotopes, precursor, product_sch)
            fn_push(data, idp, chemical_vec, mass_vec, abundance_vec, preab_vec, id)
        end
    end
    chemical = ChainedVector(data[1])
    mass_product = ChainedVector(data[2])
    abundance_pair = ChainedVector(data[3])
    id_pair = ChainedVector([[i for _ in eachindex(x)] for (i, x) in enumerate(data[1])])
    abtype = abtyped(abtype)
    if isempty(abundance_pair) 
        if iter
            (; ID = id_pair, Chemical = ChemicalTransition[], Mass = mass_product, Abundance = abundance_pair, Preab = ChainedVector(data[4]))
        else
            (; ID = id_pair, Chemical = ChemicalTransition[], Mass = mass_product, Abundance = abundance_pair)
        end
    else
        if iter 
            (; ID = id_pair, Chemical = chemical, Mass = mass_product, Abundance = abundance_pair, Preab = ChainedVector(data[4]))
        else
            (; ID = id_pair, Chemical = chemical, Mass = mass_product, Abundance = abundance_pair) 
        end
    end
end

function isotopologues_elements_single_ms2(precursor_dictionary::Dict, first_abundance, element_product, isotopes, nisotopes, iid, msfix, threshold, precise, abs_charge)
    first_element = element_product[iid]
    first_mass = (deltammi(isotopes, element_product) + msfix) / abs_charge
    if first_abundance < threshold
        element_vec = typeof(first_element)[]
        mass_vec = typeof(first_mass)[]
        abundance_vec = typeof(first_abundance)[]
    else
        element_vec = [first_element]
        mass_vec = [first_mass]
        abundance_vec = [first_abundance]
    end
    rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_product, precursor_dictionary, isotopes, nisotopes, iid, 2, 1, 1, first_mass, first_abundance, threshold, (true, true), precise, abs_charge)
    element_vec, mass_vec, abundance_vec, nothing
end

function isotopologues_elements_single_ms2_iter(precursor_dictionary::Dict, first_abundance, element_product, isotopes, nisotopes, iid, msfix, threshold, precise, abs_charge)
    first_element = element_product[iid]
    first_mass = (deltammi(isotopes, element_product) + msfix) / abs_charge
    first_preab = isotopologue_inverse_combination(precise, element_product, isotopes)
    if first_abundance < threshold
        element_vec = typeof(first_element)[]
        mass_vec = typeof(first_mass)[]
        abundance_vec = typeof(first_abundance)[]
        preab_vec = typeof(first_preab)[]
    else
        element_vec = [first_element]
        mass_vec = [first_mass]
        abundance_vec = [first_abundance]
        preab_vec = [first_preab]
    end
    rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_product, precursor_dictionary, isotopes, nisotopes, iid, 2, 1, 1, first_mass, first_abundance, first_preab, threshold, (true, true), precise, abs_charge)
    element_vec, mass_vec, abundance_vec, preab_vec
end

function product_id_sort(it1, idp, element_vec, mass_vec, isotopes, precursor_isotopes, product_isotopes, precursor, product_sch) 
    id = sortperm(mass_vec) 
    [vcat(precursor[idp], isotopomerize(product_sch, ElementsVector(product_isotopes, element_vec[i]))) for i in id], id 
end
product_id(it1, idp, element_vec, mass_vec, isotopes, precursor_isotopes, product_isotopes, precursor, product_sch) = [vcat(precursor[idp], isotopomerize(product_sch, ElementsVector(product_isotopes, element))) for element in element_vec], nothing

function loss_id_sort(it1, idp, element_vec, mass_vec, isotopes, precursor_isotopes, product_isotopes, precursor, product_sch)  
    id = sortperm(mass_vec)
    element_vec = number_loss_elements(isotopes, element_vec, precursor_isotopes, precursor_loss(precursor_isotopes, it1.Element[idp]), id)
    [vcat(precursor[idp], isotopomerize(product_sch, ElementsVector(product_isotopes, element))) for element in element_vec], id
end

function loss_id(it1, idp, element_vec, mass_vec, isotopes, precursor_isotopes, product_isotopes, precursor, product_sch)
    element_vec = number_loss_elements(isotopes, element_vec, precursor_isotopes, precursor_loss(precursor_isotopes, it1.Element[idp]), eachindex(mass_vec))
    [vcat(precursor[idp], isotopomerize(product_sch, ElementsVector(product_isotopes, element))) for element in element_vec], nothing
end

function precursor_loss(precursor_isotopes, elements)
    precursor_vec = zeros(Int, length(precursor_isotopes))
    for (k, v) in elements
        i = findfirst(==(k), precursor_isotopes)
        isnothing(i) && continue
        precursor_vec[i] = -v 
    end
    precursor_vec
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

function push_data_sort!(data, idp, chemical_vec, mass_vec, abundance_vec, preab_vec, id)
    data[1][idp] = chemical_vec
    data[2][idp] = mass_vec[id]
    data[3][idp] = abundance_vec[id]
    data
end

function push_data!(data, idp, chemical_vec, mass_vec, abundance_vec, preab_vec, id)
    data[1][idp] = chemical_vec
    data[2][idp] = mass_vec
    data[3][idp] = abundance_vec
    data
end

function push_data_iter_sort!(data, idp, chemical_vec, mass_vec, abundance_vec, preab_vec, id)
    data[1][idp] = chemical_vec
    data[2][idp] = mass_vec[id]
    data[3][idp] = abundance_vec[id]
    data[4][idp] = preab_vec[id]
    data
end

function push_data_iter!(data, idp, chemical_vec, mass_vec, abundance_vec, preab_vec, id)
    data[1][idp] = chemical_vec
    data[2][idp] = mass_vec
    data[3][idp] = abundance_vec
    data[4][idp] = preab_vec
    data
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
        precise,
        abs_charge
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
        next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, prev_mass, prev_proportion, threshold, current, precise, abs_charge)
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
            new_proportion = update_proportion(precise, proportion, rel, ris)
            new_proportion = update_proportion(precise, new_proportion, pis, pel)
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
            mass += element_mass_delta(i, e) / abs_charge
            if !any(current)
                backward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                push!(element_vec, element_vp[iid])
                push!(abundance_vec, proportion)
                push!(mass_vec, mass) 
                backward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, threshold, next, precise, abs_charge)
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
                backward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, threshold, next, precise, abs_charge)
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
            new_proportion = update_proportion(precise, proportion, pel, pis)
            new_proportion = update_proportion(precise, new_proportion, ris, rel)
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
            mass += element_mass_delta(e, i) / abs_charge
            if !any(current)
                forward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                push!(element_vec, element_vp[iid])
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, threshold, next, precise, abs_charge)
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
                forward_next_proportion = rec_exchangeisotopes!(element_vec, mass_vec, abundance_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, threshold, next, precise, abs_charge)
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
        precise,
        abs_charge
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
        next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, prev_mass, prev_proportion, prev_preab, threshold, current, precise, abs_charge)
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
            new_proportion = update_proportion(precise, proportion, rel, ris)
            new_proportion = update_proportion(precise, new_proportion, pis, pel)
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
            preab = update_inverse_proportion(precise, preab, pis, pel)
            pis -= 1
            pel += 1
            ris += 1 
            rel -= 1
            proportion = new_proportion
            mass += element_mass_delta(i, e) / abs_charge
            if !any(current)
                backward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass) 
                backward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, preab, threshold, next, precise, abs_charge)
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
                backward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, preab, threshold, next, precise, abs_charge)
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
            new_proportion = update_proportion(precise, proportion, pel, pis)
            new_proportion = update_proportion(precise, new_proportion, ris, rel)
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
            preab = update_inverse_proportion(precise, preab, pel, pis)
            pis += 1
            pel -= 1
            ris -= 1 
            rel += 1
            proportion = new_proportion
            mass += element_mass_delta(e, i) / abs_charge
            if !any(current)
                forward_max_proportion = proportion
            elseif proportion >= threshold && iter
                element_vp[isotope_position] = pis
                element_vp[element_position] = pel
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, proportion)
                push!(mass_vec, mass)
                forward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, preab, threshold, next, precise, abs_charge)
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
                forward_next_proportion = rec_exchangeisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, element_precursor_dictionary, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, proportion, preab, threshold, next, precise, abs_charge)
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