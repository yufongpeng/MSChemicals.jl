# ==========================================================================================================================
# Mid level MSn
function isotopologues_elements_msn(precise::Val, element_vp_vec, msfix_vec, abundance, abtype, threshold) 
    max_vp_vec = map(x -> maximal_abundance_elements_composition(precise, x), element_vp_vec)
    # max_proportion_vec = map(x -> isotopicabundance(precise, x), max_dictionary_vec) 
    max_proportion = [prod(first(x)) for x in max_vp_vec]
    total, abundance_cutoff, proportion_cutoff = abundance_threshold_vec(abtype, abundance, threshold, max_proportion, element_vp_vec) 
    tbls = map(isotopologues_elements_msx, element_vp_vec, msfix_vec, max_vp_vec, [proportion_cutoff for _ in eachindex(element_vp_vec)], [precise for _ in eachindex(element_vp_vec)])
    max_proportion[begin] *= total
    els = Vector{ElementsVector}[]
    abv = typeof(return_abundance(precise, big(0.0)))[]
    mass = Vector{float(Int)}[]
    combinemsx!(els, mass, abv, tbls, abundance_cutoff, prod(max_proportion), [1 for _ in eachindex(tbls)], 1)
    # rec_vec_ab!(els, abv, mass, tbls, abundance_cutoff, prod(max_proportion), [1 for _ in eachindex(tbls)], 1)
    idm = sortperm(mass)
    abv = abv[idm]
    if dopostnormalize(abtype)
        abv = normalize_abundance(abv, abundance, abtype)
    end
    (; Element = parallel_gain_elements!(els[idm]), Mass = mass[idm], Abundance = abv)
end

function isotopologues_elements_msx(element_vp, msfix, max_vp, proportion_cutoff, precise)
    isempty(element_vp) && return (; Element = [ElementsVector(String[], Int[])], Mass = [mmi(element_vp) + msfix], Abundance = [float(1)])
    isotopes, els, mass, abv, _ = isotopologues_elements_single(precise, element_vp, msfix, 1.0, Max(), rcrit(proportion_cutoff), max_vp..., false)
    id = sortperm(abv; rev = true)
    (; Element = [ElementsVector(isotopes, els[i]) for i in id], Mass = mass[id], Abundance = abv[id]) 
end

# ==========================================================================================================================
# Low level MSn
function combinemsx!(els, mass, abv, tbls, abundance_cutoff, maxab, id, eln)
    eln == lastindex(tbls) && return combinemsx_end!(els, mass, abv, tbls, abundance_cutoff, maxab, id)
    # maxab /= first(tbls[eln].Abundance)
    # pass = false
    i = 0
    @inbounds while i < length(tbls[eln].Abundance)
        id[eln] = i + 1
        combinemsx!(els, mass, abv, tbls, abundance_cutoff, maxab * tbls[eln].Abundance[i + 1], id, eln + 1) || break
        i += 1
    end
    id[eln] = 1
    i > 0
end

function combinemsx_end!(els, mass, abv, tbls, abundance_cutoff, maxab, id)
    i = 0
    @inbounds while i < length(tbls[end].Abundance)
        ab = maxab * tbls[end].Abundance[i + 1]
        ab < abundance_cutoff && break
        i += 1
        id[end] = i
        push!(els, [tbl.Element[i] for (tbl, i) in zip(tbls, id)])
        push!(abv, ab)
        push!(mass, reverse(cumsum(tbl.Mass[i] for (tbl, i) in Iterators.reverse(zip(tbls, id)))))
    end
    id[end] = 1
    i > 0
end