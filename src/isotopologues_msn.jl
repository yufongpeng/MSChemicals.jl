# ==========================================================================================================================
# Mid level MSn
function isotopologues_elements_msn(precise::Val, element_vp_vec, msfix_vec, abundance, abtype, threshold, sort) 
    max_vp_vec = map(x -> maximal_abundance_elements_composition(precise, x), element_vp_vec)
    max_proportion = [prod(first(x)) for x in max_vp_vec]
    total, abundance_cutoff, proportion_cutoff = abundance_threshold_vec(abtype, abundance, threshold, max_proportion, element_vp_vec) 
    tbls = map(isotopologues_elements_msx, element_vp_vec, msfix_vec, max_vp_vec, [proportion_cutoff for _ in eachindex(element_vp_vec)], [precise for _ in eachindex(element_vp_vec)], [1 for _ in eachindex(element_vp_vec)])
    max_proportion[begin] *= total
    els, mass, abv = combinemsx(precise, tbls, abundance_cutoff, prod(max_proportion))
    if sort
        idm = sortperm(mass)
        abv = abv[idm]
        els = els[idm]
        mass = mass[idm]
    end
    if dopostnormalize(abtype)
        abv = normalize_abundance(abv, abundance, abtype)
    end
    (; Element = parallel_gain_elements!(els), Mass = mass, Abundance = abv)
end

function isotopologues_elements_msx(element_vp, msfix, max_vp, proportion_cutoff, precise, abs_charge)
    isempty(element_vp) && return (; Element = [ElementsVector(String[], Int[])], Mass = [mmi(element_vp) + msfix], Abundance = [float(1)])
    isotopes, els, mass, abv, _ = isotopologues_elements_single(precise, element_vp, msfix, 1.0, Max(), rcrit(proportion_cutoff), max_vp..., false, abs_charge)
    id = sortperm(abv; rev = true)
    (; Element = [ElementsVector(isotopes, els[i]) for i in id], Mass = mass[id], Abundance = abv[id]) 
end

# ==========================================================================================================================
# Low level MSn
function combinemsx(precise, tbls, abundance_cutoff, maxab)
    els = Vector{ElementsVector}[]
    mass = Vector{float(Int)}[]
    abv = typeof(return_abundance(precise, big(0.0)))[]
    el = [tbl.Element[begin] for tbl in tbls]
    ms = reverse(cumsum(tbl.Mass[begin] for tbl in Iterators.reverse(tbls)))
    @inbounds for i in eachindex(tbls[end].Abundance)
        push!(els, copy(el))
        els[end][end] = tbls[end].Element[i]
        push!(abv, maxab * tbls[end].Abundance[i])
        push!(mass, copy(ms))
        mass[end][end] = tbls[end].Mass[i]
    end
    ci = length(tbls[end].Abundance)
    @inbounds for eln in (lastindex(tbls) - 1):-1:1
        ni = ci
        for i in eachindex(tbls[eln].Abundance)
            if i == 1 
                for j in 1:ci
                    mass[j][eln] = tbls[eln].Mass[i] + mass[j][eln + 1]
                end
                continue
            end
            for j in 1:ci
                if abv[j] <= 0 
                    continue
                else
                    ab = abv[j] * tbls[eln].Abundance[i]
                end
                if ab < abundance_cutoff 
                    continue
                else
                    ni += 1
                    push!(els, copy(els[j]))
                    els[ni][eln] = tbls[eln].Element[i]
                    push!(abv, ab)
                    push!(mass, copy(mass[j]))
                    mass[ni][eln] = tbls[eln].Mass[i] + mass[j][eln + 1]
                end
            end
        end
        ci = ni
    end
    els, mass, abv
end