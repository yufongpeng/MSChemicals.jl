# ==========================================================================================================================
# Mid level MS1
isotopologues_elements(precise::Val, x::AbstractString, abundance, abtype, threshold; normalize = true) = 
    isotopologues_elements(precise, chemicalelements(x), abundance, abtype, threshold; normalize)
isotopologues_elements(precise::Val, input_element::Vector, abundance, abtype, threshold; normalize = true) = 
    isotopologues_elements(precise, get_element_dictionary_fixmass(input_element)..., abundance, abtype, threshold; normalize)
function isotopologues_elements(precise::Val, element_dictionary::Dict, msfix, abundance, abtype, threshold; normalize = true)
    abtype = abtyped(abtype)
    isempty(element_dictionary) && return (; Element = [Pair{String, Int}[]], Mass = [mmi(element_dictionary)], Abundance = [abundance], Preab = [one(abundance)]) 
    element_vp = collect(element_dictionary)
    first_proportion = isotopicabundance(precise, element_vp)
    max_proportion, max_vp = maximal_abundance_composition(precise, element_vp)
    # max_proportion = isotopicabundance(precise, max_dictionary)
    if max_proportion < first_proportion || max_proportion == 0
        if precise == Val(false)
            first_proportion = return_abundance(precise, isotopicabundance(Val(true), element_vp))
            max_proportion, max_vp = maximal_abundance_composition(Val(true), element_vp)
            max_proportion = return_abundance(precise, max_proportion)
        else max_proportion < first_proportion 
            throw(ArgumentError("The input chemical is too large to estimete isotopic abundance correctly."))
        end
    end
    total, abundance_cutoff = abundance_threshold(abtype, abundance, threshold, first_proportion, max_proportion)
    if abtype == Input() && first_proportion / max_proportion < min_propotion() 
        throw(ArgumentError("Isotopic abundance of input chemical is too small; try use `abtype` other than `Input()` or larger threshold`"))
    end
    nisotopes = get_nisotopes(element_vp)
    isotopes = get_isotopes(element_vp)
    iid = findall(!ismajor, isotopes)
    element_chemical = [max_vp[iid]]
    mass_chemical = [nmmi(isotopes, max_vp) + msfix]
    abundance_chemical = [total * max_proportion]
    rec_addminusisotopes!(element_chemical, mass_chemical, abundance_chemical, max_vp, isotopes, nisotopes, iid, 2, 1, 1, first(mass_chemical), first(abundance_chemical), abundance_cutoff, (true, true), precise)
    # abundance_chemical_normalize = normalize_abundance(abundance_chemical, abundance, preabtype(abtype), [Max(), Input(), Total()])
    id = sortperm(mass_chemical)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(abundance_chemical)))
    filter!(x -> >=(abundance_chemical[x], abundance_cutoff), id)
    if normalize && dopostnormalize(abtype)
        abundance_chemical = normalize_abundance(abundance_chemical[id], abundance, abtype, [Max(), Input(), List(), Total()])
    else
        abundance_chemical = abundance_chemical[id]
    end
    isotopes = isotopes[iid]
    # (; Element = element_chemical[id], Mass = mass_chemical[id], Abundance = abundance_chemical) 
    (; Element = [ElementsVector(isotopes, element_chemical[i]) for i in id], Mass = mass_chemical[id], Abundance = abundance_chemical) 
end

isotopologues_elements_iter(precise::Val, x::AbstractString, abundance, abtype, threshold; normalize = true) = 
    isotopologues_elements_iter(precise, precise::Val, chemicalelements(x), abundance, abtype, threshold; normalize)
isotopologues_elements_iter(precise::Val, input_element::Vector, abundance, abtype, threshold; normalize = true) = 
    isotopologues_elements_iter(precise, get_element_dictionary_fixmass(input_element)..., abundance, abtype, threshold; normalize)
function isotopologues_elements_iter(precise::Val, element_dictionary::Dict, msfix, abundance, abtype, threshold; normalize = true)
    abtype = abtyped(abtype)
    isempty(element_dictionary) && return (; Element = [Pair{String, Int}[]], Mass = [mmi(element_dictionary)], Abundance = [abundance], Preab = [one(abundance)]) 
    element_vp = collect(element_dictionary)
    first_proportion = isotopicabundance(precise, element_vp)
    max_proportion, max_vp = maximal_abundance_composition(precise, element_vp)
    # max_proportion = isotopicabundance(precise, max_dictionary)
    if max_proportion < first_proportion || max_proportion == 0
        if precise == Val(false)
            first_proportion = return_abundance(precise, isotopicabundance(Val(true), element_vp))
            max_proportion, max_vp = maximal_abundance_composition(Val(true), element_vp)
            max_proportion = return_abundance(precise, max_proportion)
        else max_proportion < first_proportion 
            throw(ArgumentError("The input chemical is too large to estimete isotopic abundance correctly."))
        end
    end
    total, abundance_cutoff = abundance_threshold(abtype, abundance, threshold, first_proportion, max_proportion)
    if abtype == Input() && first_proportion / max_proportion < min_propotion() 
        throw(ArgumentError("Isotopic abundance of input chemical is too small; try use `abtype` other than `Input()` or larger threshold`"))
    end
    nisotopes = get_nisotopes(element_vp)
    isotopes = get_isotopes(element_vp)
    iid = findall(!ismajor, isotopes)
    element_chemical = [max_vp[iid]]
    mass_chemical = [nmmi(isotopes, max_vp) + msfix]
    abundance_chemical = [total * max_proportion]
    preab_chemical = [isotopologue_inverse_combination(precise, max_vp, isotopes)]
    # rec_addminusisotopes_iter!(element_chemical, mass_chemical, abundance_chemical, preab_chemical, max_dictionary, element_isotope_pairs(element_dictionary; sort = false), 1, first(mass_chemical), first(abundance_chemical), first(preab_chemical), abundance_cutoff, true, true, precise)
    rec_addminusisotopes_iter!(element_chemical, mass_chemical, abundance_chemical, preab_chemical, max_vp, isotopes, nisotopes, iid, 2, 1, 1, first(mass_chemical), first(abundance_chemical), first(preab_chemical), abundance_cutoff, (true, true), precise)
    # abundance_chemical_normalize = normalize_abundance(abundance_chemical, abundance, preabtype(abtype), [Max(), Input(), Total()])
    id = sortperm(mass_chemical)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(abundance_chemical)))
    filter!(x -> >=(abundance_chemical[x], abundance_cutoff), id)
    if normalize && dopostnormalize(abtype)
        abundance_chemical = normalize_abundance(abundance_chemical[id], abundance, abtype, [Max(), Input(), List(), Total()])
    else
        abundance_chemical = abundance_chemical[id]
    end
    isotopes = isotopes[iid]
    # (; Element = element_chemical[id], Mass = mass_chemical[id], Abundance = abundance_chemical, Preab = preab_chemical[id]) 
    (; Element = [ElementsVector(isotopes, element_chemical[i]) for i in id], Mass = mass_chemical[id], Abundance = abundance_chemical, Preab = preab_chemical[id]) 
end

# ==========================================================================================================================
# Low level MS1
function rec_addminusisotopes!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        element_vp::Vector, 
        isotopes::Vector,
        nisotopes::Vector, 
        iid::Vector, 
        isotope_position::Int, 
        element_position::Int,
        n_position::Int,
        prev_mass,
        prev_abundance, 
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
        n_position > lastindex(nisotopes) && return prev_abundance
        element_position += 1
        isotope_position += 1
    end
    isotope_position > lastindex(isotopes) && return prev_abundance
    max_abundance = prev_abundance
    next_abundance = prev_abundance
    iter = isotope_position < lastindex(isotopes)
    iter_e = iter && (isotope_position < element_position + nisotopes[n_position] - 1)
    if iter && (prev_abundance >= threshold || all(current))
        next_abundance = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_vp, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, prev_mass, prev_abundance, threshold, current, precise)
        max_abundance = max(max_abundance, next_abundance)
    end
    ne = element_vp[element_position]
    ni = element_vp[isotope_position]
    e = isotopes[element_position]
    i = isotopes[isotope_position]
    if first(current)
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        backward_max_abundance = max(prev_abundance, next_abundance)
        while cni > 0
            new_abundance = update_abundance(precise, abundance, i, e, cni, cne, 1)
            if new_abundance >= threshold
                next = (true, true)
            elseif new_abundance >= abundance
                next = (iter_e && backward_max_abundance >= threshold) ? (false, true) : (false, false)
            elseif iter_e
                backward_max_abundance >= threshold || break
                next = (false, true)
            else
                break
            end
            cni -= 1
            cne += 1
            mass += element_mass_delta(i, e)
            abundance = new_abundance
            if !any(next)
                backward_max_abundance = abundance
            elseif abundance >= threshold && iter
                element_vp[isotope_position] = cni
                element_vp[element_position] = cne
                push!(element_vec, element_vp[iid])
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                backward_next_abundance = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_vp, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, abundance, threshold, next, precise)
                backward_max_abundance = max(abundance, backward_next_abundance)
            elseif abundance >= threshold
                element_vp[isotope_position] = cni
                push!(element_vec, element_vp[iid])
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                backward_max_abundance = abundance
            elseif iter 
                element_vp[isotope_position] = cni
                element_vp[element_position] = cne
                backward_next_abundance = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_vp, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, abundance, threshold, next, precise)
                backward_max_abundance = max(abundance, backward_next_abundance)
            else
                backward_max_abundance = abundance
            end
            max_abundance = max(backward_max_abundance, max_abundance)
        end
        element_vp[isotope_position] = ni
        element_vp[element_position] = ne
    end
    if last(current)
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        forward_max_abundance = max(prev_abundance, next_abundance)
        while cne > 0
            new_abundance = update_abundance(precise, abundance, e, i, cne, cni, 1)
            if new_abundance >= threshold
                next = (true, true)
            elseif new_abundance >= abundance
                # newi = last(element_isotope_pair[isotope_position + 1])
                next = (iter_e && forward_max_abundance >= threshold) ? (false, true) : (false, false)
            elseif iter_e
                # newi = last(element_isotope_pair[isotope_position + 1])
                # get(element_vp, newi, 0) == 0 && break 
                forward_max_abundance >= threshold || break
                next = (true, false)
            else
                break
            end
            cni += 1
            cne -= 1
            mass += element_mass_delta(e, i)
            abundance = new_abundance
            if !any(next)
                forward_max_abundance = abundance
            elseif abundance >= threshold && iter
                element_vp[isotope_position] = cni
                element_vp[element_position] = cne
                push!(element_vec, element_vp[iid])
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                forward_next_abundance = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_vp, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, abundance, threshold, next, precise)
                forward_max_abundance = max(abundance, forward_next_abundance)
            elseif abundance >= threshold
                element_vp[isotope_position] = cni
                push!(element_vec, element_vp[iid])
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                forward_max_abundance = abundance
            elseif iter
                element_vp[isotope_position] = cni
                element_vp[element_position] = cne
                forward_next_abundance = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_vp, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, abundance, threshold, next, precise)
                forward_max_abundance = max(abundance, forward_next_abundance)
            else
                forward_max_abundance = abundance
            end
            max_abundance = max(forward_max_abundance, max_abundance)
        end
        element_vp[isotope_position] = ni
        element_vp[element_position] = ne
    end
    max_abundance
end

function rec_addminusisotopes_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_vp::Vector, 
        isotopes::Vector,
        nisotopes::Vector, 
        iid::Vector, 
        isotope_position::Int, 
        element_position::Int,
        n_position::Int,
        prev_mass,
        prev_abundance, 
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
        n_position > lastindex(nisotopes) && return prev_abundance
        element_position += 1
        isotope_position += 1
    end
    isotope_position > lastindex(isotopes) && return prev_abundance
    max_abundance = prev_abundance
    next_abundance = prev_abundance
    iter = isotope_position < lastindex(isotopes)
    iter_e = iter && (isotope_position < element_position + nisotopes[n_position] - 1)
    if iter && (prev_abundance >= threshold || all(current))
        next_abundance = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, prev_mass, prev_abundance, prev_preab, threshold, current, precise)
        max_abundance = max(max_abundance, next_abundance)
    end
    ne = element_vp[element_position]
    ni = element_vp[isotope_position]
    e = isotopes[element_position]
    i = isotopes[isotope_position]
    if first(current)
        preab = prev_preab
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        backward_max_abundance = max(prev_abundance, next_abundance)
        while cni > 0
            new_abundance = update_abundance(precise, abundance, i, e, cni, cne, 1)
            if new_abundance >= threshold
                next = (true, true)
            elseif new_abundance >= abundance
                next = (iter_e && backward_max_abundance >= threshold) ? (false, true) : (false, false)
            elseif iter_e
                backward_max_abundance >= threshold || break
                next = (false, true)
            else
                break
            end
            preab = update_inverse_proportion(precise, preab, cni, cne, 1)
            cni -= 1
            cne += 1
            mass += element_mass_delta(i, e)
            abundance = new_abundance
            if !any(next)
                backward_max_abundance = abundance
            elseif abundance >= threshold && iter
                element_vp[isotope_position] = cni
                element_vp[element_position] = cne
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                backward_next_abundance = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, abundance, preab, threshold, next, precise)
                backward_max_abundance = max(abundance, backward_next_abundance)
            elseif abundance >= threshold
                element_vp[isotope_position] = cni
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                backward_max_abundance = abundance
            elseif iter 
                element_vp[isotope_position] = cni
                element_vp[element_position] = cne
                backward_next_abundance = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, abundance, preab, threshold, next, precise)
                backward_max_abundance = max(abundance, backward_next_abundance)
            else
                backward_max_abundance = abundance
            end
            max_abundance = max(max_abundance, backward_max_abundance)
        end
        element_vp[isotope_position] = ni
        element_vp[element_position] = ne
    end
    if last(current)
        preab = prev_preab
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        iter_c = false
        forward_max_abundance = max(prev_abundance, next_abundance)
        while cne > 0
            new_abundance = update_abundance(precise, abundance, e, i, cne, cni, 1)
            if new_abundance >= threshold
                next = (true, true)
            elseif new_abundance >= abundance
                # newi = last(element_isotope_pair[isotope_position + 1])
                next = (iter_e && forward_max_abundance >= threshold) ? (false, true) : (false, false)
            elseif iter_e
                # newi = last(element_isotope_pair[isotope_position + 1])
                # get(element_vp, newi, 0) == 0 && break 
                forward_max_abundance >= threshold || break
                next = (true, false)
            else
                break
            end
            preab = update_inverse_proportion(precise, preab, cne, cni, 1)
            cni += 1
            cne -= 1
            mass += element_mass_delta(e, i)
            abundance = new_abundance
            if !any(next)
                forward_max_abundance = abundance
            elseif abundance >= threshold && iter
                element_vp[isotope_position] = cni
                element_vp[element_position] = cne
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                forward_next_abundance = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, abundance, preab, threshold, next, precise)
                forward_max_abundance = max(abundance, forward_next_abundance)
            elseif abundance >= threshold
                element_vp[isotope_position] = cni
                push!(element_vec, element_vp[iid])
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                forward_max_abundance = abundance
            elseif iter
                element_vp[isotope_position] = cni
                element_vp[element_position] = cne
                forward_next_abundance = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, isotopes, nisotopes, iid, isotope_position + 1, element_position, n_position, mass, abundance, preab, threshold, next, precise)
                forward_max_abundance = max(abundance, forward_next_abundance)
            else
                backward_max_abundance = abundance
            end
            max_abundance = max(max_abundance, forward_max_abundance)
        end
        element_vp[isotope_position] = ni
        element_vp[element_position] = ne
    end
    max_abundance
end