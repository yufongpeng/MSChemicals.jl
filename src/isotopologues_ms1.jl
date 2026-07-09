# ==========================================================================================================================
# Mid level MS1
isotopologues_elements(precise::Val, x::AbstractString, abundance, abtype, threshold; normalize = true) = 
    isotopologues_elements(precise, chemicalelements(x), abundance, abtype, threshold; normalize)
isotopologues_elements(precise::Val, input_element::Vector, abundance, abtype, threshold; normalize = true) = 
    isotopologues_elements(precise, get_element_dictinonary_fixmass(input_element)..., abundance, abtype, threshold; normalize)
function isotopologues_elements(precise::Val, element_dictionary::Dict, msfix, abundance, abtype, threshold; normalize = true)
    abtype = abtyped(abtype)
    isempty(element_dictionary) && return (; Element = [element_dictionary], Mass = [mmi(element_dictionary)], Abundance = [abundance]) 
    first_proportion = isotopicabundance(precise, element_dictionary)
    max_dictionary = maximal_elements(element_dictionary)
    max_proportion = isotopicabundance(precise, max_dictionary)
    if max_proportion < first_proportion || max_proportion == 0
        if precise == Val(false)
            first_proportion = return_abundance(precise, isotopicabundance(Val(true), element_dictionary))
            max_proportion = return_abundance(precise, isotopicabundance(Val(true), max_dictionary))
        else max_proportion < first_proportion 
            throw(ArgumentError("The input chemical is too large to estimete isotopic abundance correctly."))
        end
    end
    total, abundance_cutoff = abundance_threshold(abtype, abundance, threshold, first_proportion, max_proportion)
    if abtype == Input() && first_proportion / max_proportion < min_propotion() 
        throw(ArgumentError("Isotopic abundance of input chemical is too small; try use `abtype` other than `Input()` or larger threshold`"))
    end
    element_chemical = [get_isotope_vec(max_dictionary)]
    mass_chemical = [mmi(max_dictionary) + msfix]
    abundance_chemical = [total * max_proportion]
    rec_addminusisotopes!(element_chemical, mass_chemical, abundance_chemical, max_dictionary, element_isotope_pairs(element_dictionary; sort = false), 1, first(mass_chemical), first(abundance_chemical), abundance_cutoff, (true, true), precise)
    # abundance_chemical_normalize = normalize_abundance(abundance_chemical, abundance, preabtype(abtype), [Max(), Input(), Total()])
    id = sortperm(mass_chemical)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(abundance_chemical)))
    filter!(x -> >=(abundance_chemical[x], abundance_cutoff), id)
    if normalize && dopostnormalize(abtype)
        abundance_chemical = normalize_abundance(abundance_chemical[id], abundance, abtype, [Max(), Input(), List(), Total()])
    else
        abundance_chemical = abundance_chemical[id]
    end
    (; Element = element_chemical[id], Mass = mass_chemical[id], Abundance = abundance_chemical) 
end

isotopologues_elements_iter(precise::Val, x::AbstractString, abundance, abtype, threshold; normalize = true) = 
    isotopologues_elements_iter(precise, precise::Val, chemicalelements(x), abundance, abtype, threshold; normalize)
isotopologues_elements_iter(precise::Val, input_element::Vector, abundance, abtype, threshold; normalize = true) = 
    isotopologues_elements_iter(precise, get_element_dictinonary_fixmass(input_element)..., abundance, abtype, threshold; normalize)
function isotopologues_elements_iter(precise::Val, element_dictionary::Dict, msfix, abundance, abtype, threshold; normalize = true)
    abtype = abtyped(abtype)
    isempty(element_dictionary) && return (; Element = [element_dictionary], Mass = [mmi(element_dictionary)], Abundance = [abundance], Preab = [one(abundance)]) 
    first_proportion = isotopicabundance(precise, element_dictionary)
    max_dictionary = maximal_elements(element_dictionary)
    max_proportion = isotopicabundance(precise, max_dictionary)
    if max_proportion < first_proportion || max_proportion == 0
        if precise == Val(false)
            first_proportion = return_abundance(precise, isotopicabundance(Val(true), element_dictionary))
            max_proportion = return_abundance(precise, isotopicabundance(Val(true), max_dictionary))
        else max_proportion < first_proportion 
            throw(ArgumentError("The input chemical is too large to estimete isotopic abundance correctly."))
        end
    end
    total, abundance_cutoff = abundance_threshold(abtype, abundance, threshold, first_proportion, max_proportion)
    if abtype == Input() && first_proportion / max_proportion < min_propotion() 
        throw(ArgumentError("Isotopic abundance of input chemical is too small; try use `abtype` other than `Input()` or larger threshold`"))
    end
    element_chemical = [get_isotope_vec(max_dictionary)]
    mass_chemical = [mmi(max_dictionary) + msfix]
    abundance_chemical = [total * max_proportion]
    preab_chemical = [isotopologue_inverse_combination(precise, max_dictionary)]
    # rec_addminusisotopes_iter!(element_chemical, mass_chemical, abundance_chemical, preab_chemical, max_dictionary, element_isotope_pairs(element_dictionary; sort = false), 1, first(mass_chemical), first(abundance_chemical), first(preab_chemical), abundance_cutoff, true, true, precise)
    rec_addminusisotopes_iter!(element_chemical, mass_chemical, abundance_chemical, preab_chemical, max_dictionary, element_isotope_pairs(element_dictionary; sort = false), 1, first(mass_chemical), first(abundance_chemical), first(preab_chemical), abundance_cutoff, (true, true), precise)
    # abundance_chemical_normalize = normalize_abundance(abundance_chemical, abundance, preabtype(abtype), [Max(), Input(), Total()])
    id = sortperm(mass_chemical)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(abundance_chemical)))
    filter!(x -> >=(abundance_chemical[x], abundance_cutoff), id)
    if normalize && dopostnormalize(abtype)
        abundance_chemical = normalize_abundance(abundance_chemical[id], abundance, abtype, [Max(), Input(), List(), Total()])
    else
        abundance_chemical = abundance_chemical[id]
    end
    (; Element = element_chemical[id], Mass = mass_chemical[id], Abundance = abundance_chemical, Preab = preab_chemical[id]) 
end

# ==========================================================================================================================
# Low level MS1
function rec_addminusisotopes!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        element_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_abundance, 
        threshold,
        current,
        precise
    )
    next_state = false
    (e, i) = element_isotope_pair[isotope_position]
    iter = isotope_position < lastindex(element_isotope_pair)
    iter_e = iter && (first(element_isotope_pair[isotope_position + 1]) == e)
    if iter && (prev_abundance >= threshold || all(current))
        next_state = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_abundance, threshold, current, precise)
    end
    current_state = next_state
    ne = get(element_dictionary, e, 0)
    ni = get(element_dictionary, i, 0)
    if first(current)
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        iter_c = false
        backward_next_state = next_state
        while cni > 0
            new_abundance = update_abundance(precise, abundance, i, e, cni, cne, 1)
            if new_abundance >= threshold
                next = (true, true)
            elseif iter_e && new_abundance >= abundance
                newi = last(element_isotope_pair[isotope_position + 1])
                iter_c = update_abundance(precise, 1, e, newi, cne + 1, get(element_dictionary, newi, 0), 1) > 1 
                next = iter_c ? (false, true) : (false, false)
            elseif new_abundance >= abundance
                next = (false, false)
            elseif iter_e && iter_c && backward_next_state
                next = (false, true)
            elseif iter_e && iter_c
                break
            elseif iter_e
                newi = last(element_isotope_pair[isotope_position + 1])
                update_abundance(precise, 1, e, newi, cne + 1, get(element_dictionary, newi, 0), 1) > 1 || break
                iter_c = true
                next = (false, true)
            else
                break
            end
            cni -= 1
            cne += 1
            mass += element_mass_delta(i, e)
            abundance = new_abundance
            if !any(next)
                continue
            elseif abundance >= threshold && iter
                element_dictionary[i] = cni
                element_dictionary[e] = cne
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                current_state = true
                backward_next_state = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position + 1, mass, abundance, threshold, next, precise)
            elseif abundance >= threshold
                element_dictionary[i] = cni
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                current_state = true
            elseif iter 
                element_dictionary[i] = cni
                element_dictionary[e] = cne
                backward_next_state = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position + 1, mass, abundance, threshold, next, precise)
                current_state = current_state || backward_next_state
            end
        end
        element_dictionary[i] = ni
        element_dictionary[e] = ne
    end
    if last(current)
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        iter_c = false
        forward_next_state = next_state
        while cne > 0
            new_abundance = update_abundance(precise, abundance, e, i, cne, cni, 1)
            if new_abundance >= threshold
                next = (true, true)
            elseif iter_e && new_abundance >= abundance
                newi = last(element_isotope_pair[isotope_position + 1])
                if get(element_dictionary, newi, 0) > 0 
                    iter_c = update_abundance(precise, 1, newi, e, get(element_dictionary, newi, 0), cne - 1, 1) > 1 
                end
                next = iter_c ? (false, true) : (false, false)
            elseif new_abundance >= abundance
                next = (false, false)
            elseif iter_e && iter_c && forward_next_state
                next = (false, true)
            elseif iter_e && iter_c
                break
            elseif iter_e
                newi = last(element_isotope_pair[isotope_position + 1])
                get(element_dictionary, newi, 0) == 0 && break 
                update_abundance(precise, 1, newi, e, get(element_dictionary, newi, 0), cne - 1, 1) > 1 || break
                iter_c = true
                next = (true, false)
            else
                break
            end
            cni += 1
            cne -= 1
            mass += element_mass_delta(e, i)
            abundance = new_abundance
            if !any(next)
                continue
            elseif abundance >= threshold && iter
                element_dictionary[i] = cni
                element_dictionary[e] = cne
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                current_state = true
                forward_next_state = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position + 1, mass, abundance, threshold, next, precise)
            elseif abundance >= threshold
                element_dictionary[i] = cni
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                current_state = true
            elseif iter
                element_dictionary[i] = cni
                element_dictionary[e] = cne
                forward_next_state = rec_addminusisotopes!(element_vec, mass_vec, abundance_vec, element_dictionary, element_isotope_pair, isotope_position + 1, mass, abundance, threshold, next, precise)
                current_state = current_state || forward_next_state
            end
        end
        element_dictionary[i] = ni
        element_dictionary[e] = ne
    end
    current_state
end

function rec_addminusisotopes_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_dictionary::Dict, 
        element_isotope_pair::Vector, 
        isotope_position::Int, 
        prev_mass,
        prev_abundance, 
        prev_preab, 
        threshold, 
        current,
        precise
    )
    next_state = false
    (e, i) = element_isotope_pair[isotope_position]
    iter = isotope_position < lastindex(element_isotope_pair)
    iter_e = iter && first(element_isotope_pair[isotope_position + 1]) == e
    if iter && (prev_abundance >= threshold || all(current))
        next_state = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_dictionary, element_isotope_pair, isotope_position + 1, prev_mass, prev_abundance, prev_preab, threshold, current, precise)
    end
    current_state = next_state
    ne = get(element_dictionary, e, 0)
    ni = get(element_dictionary, i, 0)
    if first(current)
        preab = prev_preab
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        iter_c = false
        backward_next_state = next_state
        while cni > 0
            new_abundance = update_abundance(precise, abundance, i, e, cni, cne, 1)
            if new_abundance >= threshold
                next = (true, true)
            elseif iter_e && new_abundance >= abundance
                newi = last(element_isotope_pair[isotope_position + 1])
                iter_c = update_abundance(precise, 1, e, newi, cne + 1, get(element_dictionary, newi, 0), 1) > 1 
                next = iter_c ? (false, true) : (false, false)
            elseif new_abundance >= abundance
                next = (false, false)
            elseif iter_e && iter_c && backward_next_state
                next = (false, true)
            elseif iter_e && iter_c
                break
            elseif iter_e
                newi = last(element_isotope_pair[isotope_position + 1])
                update_abundance(precise, 1, e, newi, cne + 1, get(element_dictionary, newi, 0), 1) > 1 || break
                iter_c = true
                next = (false, true)
            else
                break
            end
            preab = update_inverse_proportion(precise, preab, cni, cne, 1)
            cni -= 1
            cne += 1
            mass += element_mass_delta(i, e)
            abundance = new_abundance
            if abundance >= threshold && iter
                element_dictionary[i] = cni
                element_dictionary[e] = cne
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                current_state = true
                rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_dictionary, element_isotope_pair, isotope_position + 1, mass, abundance, preab, threshold, next, precise)
            elseif abundance >= threshold
                element_dictionary[i] = cni
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                current_state = true
            elseif iter 
                element_dictionary[i] = cni
                element_dictionary[e] = cne
                backward_next_state = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_dictionary, element_isotope_pair, isotope_position + 1, mass, abundance, preab, threshold, next, precise)
                current_state = current_state || backward_next_state
            end
        end
        element_dictionary[i] = ni
        element_dictionary[e] = ne
    end
    if last(current)
        preab = prev_preab
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        iter_c = false
        forward_next_state = next_state
        while cne > 0
            new_abundance = update_abundance(precise, abundance, e, i, cne, cni, 1)
            if new_abundance >= threshold
                next = (true, true)
            elseif iter_e && new_abundance >= abundance
                newi = last(element_isotope_pair[isotope_position + 1])
                if get(element_dictionary, newi, 0) > 0 
                    iter_c = update_abundance(precise, 1, newi, e, get(element_dictionary, newi, 0), cne - 1, 1) > 1 
                end
                next = iter_c ? (false, true) : (false, false)
            elseif new_abundance >= abundance
                next = (false, false)
            elseif iter_e && iter_c && forward_next_state
                next = (false, true)
            elseif iter_e && iter_c
                break
            elseif iter_e
                newi = last(element_isotope_pair[isotope_position + 1])
                get(element_dictionary, newi, 0) == 0 && break 
                update_abundance(precise, 1, newi, e, get(element_dictionary, newi, 0), cne - 1, 1) > 1 || break
                iter_c = true
                next = (true, false)
            else
                break
            end
            preab = update_inverse_proportion(precise, preab, cne, cni, 1)
            cni += 1
            cne -= 1
            mass += element_mass_delta(e, i)
            abundance = new_abundance
            if abundance >= threshold && iter
                element_dictionary[i] = cni
                element_dictionary[e] = cne
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                current_state = true
                rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_dictionary, element_isotope_pair, isotope_position + 1, mass, abundance, preab, threshold, next, precise)
            elseif abundance >= threshold
                element_dictionary[i] = cni
                push!(element_vec, get_isotope_vec(element_dictionary))
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                current_state = true
            elseif iter
                element_dictionary[i] = cni
                element_dictionary[e] = cne
                forward_next_state = rec_addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_dictionary, element_isotope_pair, isotope_position + 1, mass, abundance, preab, threshold, next, precise)
                current_state = current_state || forward_next_state
            end
        end
        element_dictionary[i] = ni
        element_dictionary[e] = ne
    end
    current_state
end