# ==========================================================================================================================
# Mid level MS1
isotopologues_elements_ms1(precise::Val, input_chemical, x::AbstractString, abundance, abtype, threshold, iter, sort; net_charge = charge(input_chemical)) = 
    isotopologues_elements(precise, input_chemical,chemicalelements(x), abundance, abtype, threshold, iter, sort; net_charge)
isotopologues_elements_ms1(precise::Val, input_chemical, input_element::Vector, abundance, abtype, threshold, iter, sort; net_charge = charge(input_chemical)) = 
    isotopologues_elements_ms1(precise, input_chemical, get_element_dictionary_fixmass(input_element)..., abundance, abtype, threshold, iter, sort; net_charge)
function isotopologues_elements_ms1(precise::Val, input_chemical, element_dictionary::Dict, msfix, abundance, abtype, threshold, iter, sort; net_charge = charge(input_chemical))
    abs_charge = max(1, abs(net_charge))
    msfix -= net_charge * ME
    msfix /= abs_charge
    if isempty(element_dictionary) 
        chemical = [isotopomerize(input_chemical, ElementsVector(String[], Int[]))]
        mass = [msfix]
        abv = typeof(return_abundance(precise, big(0.0)))[abundance]
        if iter 
            preab = [one(first(abv))]
        end
    else
        element_vp = collect(element_dictionary)
        abtype = abtyped(abtype)
        max_proportion, max_vp = maximal_abundance_elements_composition_check(precise, element_vp, abtype)
        isotopes, els, mass, abv, preab = isotopologues_elements_single(precise, element_vp, msfix, abundance, abtype, threshold, max_proportion, max_vp, iter, abs_charge)
        if sort 
            idm = sortperm(mass)
            abv = abv[idm]
            chemical = [isotopomerize(input_chemical, ElementsVector(isotopes, els[i])) for i in idm]
            mass = mass[idm]
            if iter 
                preab = preab[idm]
            end
        else
            chemical = [isotopomerize(input_chemical, ElementsVector(isotopes, el)) for el in els]
        end 
        if dopostnormalize(abtype)
            abv = normalize_abundance(abv, abundance, abtype)
        end
    end
    if iter
        (; Chemical = chemical, Mass = mass, Abundance = abv, Preab = preab)
    else
        (; Chemical = chemical, Mass = mass, Abundance = abv)
    end
end

function isotopologues_elements_single(precise::Val, element_vp::Vector, msfix, abundance, abtype, threshold, max_proportion, max_vp::Vector, iter, abs_charge)
    total, abundance_cutoff, proportion_cutoff = abundance_threshold_vec(abtype, abundance, threshold, max_proportion, element_vp)
    max_proportion[begin] *= total
    els = Vector{Int}[]
    abv = typeof(return_abundance(precise, big(0.0)))[]
    mass = float(Int)[]
    if iter 
        tbls = map(subisotopologues_iter, element_vp, [1.0 for _ in eachindex(element_vp)], max_vp, [proportion_cutoff for _ in eachindex(element_vp)], [precise for _ in eachindex(element_vp)], [abs_charge for _ in eachindex(element_vp)])
        preab = float(Int)[]
        combinesubisotopologues_iter!(els, mass, abv, preab, tbls, abundance_cutoff, Int[],  msfix, prod(max_proportion), 1.0, 1)
    else
        tbls = map(subisotopologues, element_vp, [1.0 for _ in eachindex(element_vp)], max_vp, [proportion_cutoff for _ in eachindex(element_vp)], [precise for _ in eachindex(element_vp)], [abs_charge for _ in eachindex(element_vp)])
        preab = nothing
        combinesubisotopologues!(els, mass, abv, tbls, abundance_cutoff, Int[], msfix, prod(max_proportion), 1)
    end
    vcat((tbl.Isotope for tbl in tbls)...), els, mass, abv, preab
end

function subisotopologues(element_vp, max_proportion, max_vp, proportioon_cutoff, precise, abs_charge)
    isempty(element_vp) && return (; Isotope = String[], Number = Int[empty(max_vp)], Mass = [mmi(element_vp)], Abundance = [float(1)])
    isotopes = elements_isotopes()[first(element_vp)]
    if length(isotopes) < 2 
        return (; Isotope = empty(isotopes), Number = [empty(max_vp)], Mass = [nmmi(isotopes, max_vp) / abs_charge], Abundance = [max_proportion])
    end
    element_chemical = [max_vp[begin + 1:end]]
    mass_chemical = [nmmi(isotopes, max_vp) / abs_charge]
    abundance_chemical = [max_proportion]
    addminusisotopes!(element_chemical, mass_chemical, abundance_chemical, max_vp, isotopes, 2, first(mass_chemical), first(abundance_chemical), proportioon_cutoff, (true, true), precise, abs_charge)
    id = sortperm(abundance_chemical; rev = true)
    (; Isotope = isotopes[begin + 1:end], Number = element_chemical[id], Mass = mass_chemical[id], Abundance = abundance_chemical[id]) 
end

function subisotopologues_iter(element_vp, max_proportion, max_vp, proportioon_cutoff, precise, abs_charge)
    isempty(element_vp) && return (; Isotope = String[], Number = Int[empty(max_vp)], Mass = [mmi(element_vp)], Abundance = [float(1)], Preab = [float(1)])
    isotopes = elements_isotopes()[first(element_vp)]
    if length(isotopes) < 2 
        return (; Isotope = empty(isotopes), Number = [empty(max_vp)], Mass = [nmmi(isotopes, max_vp) / abs_charge], Abundance = [max_proportion], Preab = [isotopologue_inverse_combination(precise, max_vp, isotopes)])
    end
    element_chemical = [max_vp[begin + 1:end]]
    mass_chemical = [nmmi(isotopes, max_vp) / abs_charge]
    abundance_chemical = [max_proportion]
    preab_chemical = [isotopologue_inverse_combination(precise, max_vp, isotopes)]
    addminusisotopes_iter!(element_chemical, mass_chemical, abundance_chemical, preab_chemical, max_vp, isotopes, 2, first(mass_chemical), first(abundance_chemical), first(preab_chemical), proportioon_cutoff, (true, true), precise, abs_charge)
    id = sortperm(abundance_chemical; rev = true)
    (; Isotope = isotopes[begin + 1:end], Number = element_chemical[id], Mass = mass_chemical[id], Abundance = abundance_chemical[id], Preab = preab_chemical[id]) 
end

# ==========================================================================================================================
# Low level MS1
function addminusisotopes!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        element_vp::Vector, 
        isotopes::Vector,
        isotope_position::Int, 
        prev_mass,
        prev_abundance, 
        threshold,
        current,
        precise,
        abs_charge
    )
    isotope_position > lastindex(isotopes) && return prev_abundance
    iter = isotope_position < lastindex(isotopes)
    max_abundance = prev_abundance
    next_abundance = prev_abundance
    if (iter && (prev_abundance >= threshold || all(current)))
        next_abundance = addminusisotopes!(element_vec, mass_vec, abundance_vec, element_vp, isotopes, isotope_position + 1, prev_mass, prev_abundance, threshold, current, precise, abs_charge)
        max_abundance = max(max_abundance, next_abundance)
    end
    ne = first(element_vp)
    ni = element_vp[isotope_position]
    e = first(isotopes)
    i = isotopes[isotope_position]
    if first(current)
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        backward_max_abundance = max(prev_abundance, next_abundance)
        while cni > 0
            new_abundance = update_abundance(precise, abundance, i, e, cni, cne)
            if new_abundance >= threshold
                next = (true, true)
            elseif new_abundance >= abundance
                next = (iter && backward_max_abundance >= threshold) ? (false, true) : (false, false)
            elseif iter
                backward_max_abundance >= threshold || break
                next = (false, true)
            else
                break
            end
            cni -= 1
            cne += 1
            mass += element_mass_delta(i, e) / abs_charge
            abundance = new_abundance
            if !any(next)
                backward_max_abundance = abundance
            elseif abundance >= threshold && iter
                element_vp[isotope_position] = cni
                element_vp[begin] = cne
                push!(element_vec, element_vp[begin + 1:end])
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                backward_next_abundance = addminusisotopes!(element_vec, mass_vec, abundance_vec, element_vp, isotopes, isotope_position + 1, mass, abundance, threshold, next, precise, abs_charge)
                backward_max_abundance = max(abundance, backward_next_abundance)
            elseif abundance >= threshold
                element_vp[isotope_position] = cni
                push!(element_vec, element_vp[begin + 1:end])
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                backward_max_abundance = abundance
            elseif iter 
                element_vp[isotope_position] = cni
                element_vp[begin] = cne
                backward_next_abundance = addminusisotopes!(element_vec, mass_vec, abundance_vec, element_vp, isotopes, isotope_position + 1, mass, abundance, threshold, next, precise, abs_charge)
                backward_max_abundance = max(abundance, backward_next_abundance)
            else
                backward_max_abundance = abundance
            end
            max_abundance = max(backward_max_abundance, max_abundance)
        end
        element_vp[isotope_position] = ni
        element_vp[begin] = ne
    end
    if last(current)
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        forward_max_abundance = max(prev_abundance, next_abundance)
        while cne > 0
            new_abundance = update_abundance(precise, abundance, e, i, cne, cni)
            if new_abundance >= threshold
                next = (true, true)
            elseif new_abundance >= abundance
                next = (iter && forward_max_abundance >= threshold) ? (false, true) : (false, false)
            elseif iter
                forward_max_abundance >= threshold || break
                next = (true, false)
            else
                break
            end
            cni += 1
            cne -= 1
            mass += element_mass_delta(e, i) / abs_charge
            abundance = new_abundance
            if !any(next)
                forward_max_abundance = abundance
            elseif abundance >= threshold && iter
                element_vp[isotope_position] = cni
                element_vp[begin] = cne
                push!(element_vec, element_vp[begin + 1:end])
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                forward_next_abundance = addminusisotopes!(element_vec, mass_vec, abundance_vec, element_vp, isotopes, isotope_position + 1, mass, abundance, threshold, next, precise, abs_charge)
                forward_max_abundance = max(abundance, forward_next_abundance)
            elseif abundance >= threshold
                element_vp[isotope_position] = cni
                push!(element_vec, element_vp[begin + 1:end])
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                forward_max_abundance = abundance
            elseif iter
                element_vp[isotope_position] = cni
                element_vp[begin] = cne
                forward_next_abundance = addminusisotopes!(element_vec, mass_vec, abundance_vec, element_vp, isotopes, isotope_position + 1, mass, abundance, threshold, next, precise, abs_charge)
                forward_max_abundance = max(abundance, forward_next_abundance)
            else
                forward_max_abundance = abundance
            end
            max_abundance = max(forward_max_abundance, max_abundance)
        end
        element_vp[isotope_position] = ni
        element_vp[begin] = ne
    end
    max_abundance
end

function addminusisotopes_iter!(
        element_vec::Vector, 
        mass_vec::Vector,
        abundance_vec::Vector, 
        preab_vec::Vector, 
        element_vp::Vector, 
        isotopes::Vector,
        isotope_position::Int, 
        prev_mass,
        prev_abundance, 
        prev_preab,
        threshold,
        current,
        precise,
        abs_charge
    )
    isotope_position > lastindex(isotopes) && return prev_abundance
    iter = isotope_position < lastindex(isotopes)
    max_abundance = prev_abundance
    next_abundance = prev_abundance
    if (iter && (prev_abundance >= threshold || all(current)))
        next_abundance = addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, isotopes, isotope_position + 1, prev_mass, prev_abundance, prev_preab, threshold, current, precise, abs_charge)
        max_abundance = max(max_abundance, next_abundance)
    end
    ne = first(element_vp)
    ni = element_vp[isotope_position]
    e = first(isotopes)
    i = isotopes[isotope_position]
    if first(current)
        preab = prev_preab
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        backward_max_abundance = max(prev_abundance, next_abundance)
        while cni > 0
            new_abundance = update_abundance(precise, abundance, i, e, cni, cne)
            if new_abundance >= threshold
                next = (true, true)
            elseif new_abundance >= abundance
                next = (iter && backward_max_abundance >= threshold) ? (false, true) : (false, false)
            elseif iter
                backward_max_abundance >= threshold || break
                next = (false, true)
            else
                break
            end
            preab = update_inverse_proportion(precise, preab, cni, cne)
            cni -= 1
            cne += 1
            mass += element_mass_delta(i, e) / abs_charge
            abundance = new_abundance
            if !any(next)
                backward_max_abundance = abundance
            elseif abundance >= threshold && iter
                element_vp[isotope_position] = cni
                element_vp[begin] = cne
                push!(element_vec, element_vp[begin + 1:end])
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                backward_next_abundance = addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, isotopes, isotope_position + 1, mass, abundance, preab, threshold, next, precise, abs_charge)
                backward_max_abundance = max(abundance, backward_next_abundance)
            elseif abundance >= threshold
                element_vp[isotope_position] = cni
                push!(element_vec, element_vp[begin + 1:end])
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                backward_max_abundance = abundance
            elseif iter 
                element_vp[isotope_position] = cni
                element_vp[begin] = cne
                backward_next_abundance = addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, isotopes, isotope_position + 1, mass, abundance, preab, threshold, next, precise, abs_charge)
                backward_max_abundance = max(abundance, backward_next_abundance)
            else
                backward_max_abundance = abundance
            end
            max_abundance = max(backward_max_abundance, max_abundance)
        end
        element_vp[isotope_position] = ni
        element_vp[begin] = ne
    end
    if last(current)
        preab = prev_preab
        abundance = prev_abundance 
        mass = prev_mass
        cni = ni 
        cne = ne
        forward_max_abundance = max(prev_abundance, next_abundance)
        while cne > 0
            new_abundance = update_abundance(precise, abundance, e, i, cne, cni)
            if new_abundance >= threshold
                next = (true, true)
            elseif new_abundance >= abundance
                next = (iter && forward_max_abundance >= threshold) ? (false, true) : (false, false)
            elseif iter
                forward_max_abundance >= threshold || break
                next = (true, false)
            else
                break
            end
            preab = update_inverse_proportion(precise, preab, cne, cni)
            cni += 1
            cne -= 1
            mass += element_mass_delta(e, i) / abs_charge
            abundance = new_abundance
            if !any(next)
                forward_max_abundance = abundance
            elseif abundance >= threshold && iter
                element_vp[isotope_position] = cni
                element_vp[begin] = cne
                push!(element_vec, element_vp[begin + 1:end])
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                forward_next_abundance = addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, isotopes, isotope_position + 1, mass, abundance, preab, threshold, next, precise, abs_charge)
                forward_max_abundance = max(abundance, forward_next_abundance)
            elseif abundance >= threshold
                element_vp[isotope_position] = cni
                push!(element_vec, element_vp[begin + 1:end])
                push!(preab_vec, preab)
                push!(abundance_vec, abundance)
                push!(mass_vec, mass)
                forward_max_abundance = abundance
            elseif iter
                element_vp[isotope_position] = cni
                element_vp[begin] = cne
                forward_next_abundance = addminusisotopes_iter!(element_vec, mass_vec, abundance_vec, preab_vec, element_vp, isotopes, isotope_position + 1, mass, abundance, preab, threshold, next, precise, abs_charge)
                forward_max_abundance = max(abundance, forward_next_abundance)
            else
                forward_max_abundance = abundance
            end
            max_abundance = max(forward_max_abundance, max_abundance)
        end
        element_vp[isotope_position] = ni
        element_vp[begin] = ne
    end
    max_abundance
end

# ==========================================================================================================================
# Low level combination
function combinesubisotopologues!(els, mass, abv, tbls, abundance_cutoff, el, ms, maxab, eln)
    eln == lastindex(tbls) && return combinesubisotopologues_end!(els, mass, abv, tbls, abundance_cutoff, el, ms, maxab)
    i = 0
    @inbounds while i < length(tbls[eln].Abundance)
        combinesubisotopologues!(els, mass, abv, tbls, abundance_cutoff, vcat(el, tbls[eln].Number[i + 1]), ms + tbls[eln].Mass[i + 1], maxab * tbls[eln].Abundance[i + 1], eln + 1) || break
        i += 1
    end
    i > 0
end

function combinesubisotopologues_end!(els, mass, abv, tbls, abundance_cutoff, el, ms, maxab)
    i = 0
    @inbounds while i < length(tbls[end].Abundance)
        ab = maxab * tbls[end].Abundance[i + 1]
        ab < abundance_cutoff && break
        i += 1
        push!(els, vcat(el, tbls[end].Number[i]))
        push!(abv, ab)
        push!(mass, ms + tbls[end].Mass[i])
    end
    i > 0
end

function combinesubisotopologues_iter!(els, mass, abv, preabv, tbls, abundance_cutoff, el, ms, maxab, preab, eln)
    eln == lastindex(tbls) && return combinesubisotopologues_end_iter!(els, mass, abv, preabv, tbls, abundance_cutoff, el, ms, maxab, preab)
    i = 0
    @inbounds while i < length(tbls[eln].Abundance)
        combinesubisotopologues_iter!(els, mass, abv, preabv, tbls, abundance_cutoff, vcat(el, tbls[eln].Number[i + 1]), ms + tbls[eln].Mass[i + 1], maxab * tbls[eln].Abundance[i + 1], preab * tbls[eln].Preab[i + 1], eln + 1) || break
        i += 1
    end
    i > 0
end

function combinesubisotopologues_end_iter!(els, mass, abv, preabv, tbls, abundance_cutoff, el, ms, maxab, preab)
    i = 0
    @inbounds while i < length(tbls[end].Abundance)
        ab = maxab * tbls[end].Abundance[i + 1]
        ab < abundance_cutoff && break
        i += 1
        push!(els, vcat(el, tbls[end].Number[i]))
        push!(abv, ab)
        push!(mass, ms + tbls[end].Mass[i])
        push!(preabv, preab * tbls[end].Preab[i])
    end
    i > 0
end