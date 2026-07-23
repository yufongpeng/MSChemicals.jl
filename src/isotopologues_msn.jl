# ==========================================================================================================================
# Mid level MSn
function isotopologues_elements_msn(precise::Val, element_vp_vec, msfix_vec, abundance, abtype, threshold) 
    max_vp_vec = map(x -> maximal_abundance_composition(precise, x), element_vp_vec)
    # max_proportion_vec = map(x -> isotopicabundance(precise, x), max_dictionary_vec) 
    total, abundance_cutoff = abundance_threshold_msn(abtype, abundance, threshold, max_vp_vec, element_vp_vec) 
    proportioon_cutoff = abundance_cutoff / total
    tbls = map(isotopologues_elements_msx, element_vp_vec, msfix_vec, max_vp_vec, [proportioon_cutoff for _ in eachindex(element_vp_vec)], [precise for _ in eachindex(element_vp_vec)])
    first(tbls).Abundance .*= total
    els = Vector{ElementsVector}[]
    abv = float(Int)[]
    mass = Vector{float(Int)}[]
    rec_vec_ab!(els, abv, mass, tbls, abundance_cutoff, prod(first(tbl.Abundance) for tbl in tbls), [1 for _ in eachindex(tbls)], 1)
    idm = sortperm(mass)
    abv = abv[idm]
    if dopostnormalize(abtype)
        abv = normalize_abundance(abv, abundance, abtype)
    end
    (; Element = parallel_gain_elements!(els[idm]), Mass = mass[idm], Abundance = abv)
end

function isotopologues_elements_msx(element_vp, msfix, max_vp, proportioon_cutoff, precise)
    isempty(element_vp) && return (; Element = [ElementsVector(String[], Int[])], Mass = [mmi(element_vp) + msfix], Abundance = [float(1)])
    nisotopes = get_nisotopes(element_vp)
    isotopes = get_isotopes(element_vp)
    iid = findall(!ismajor, isotopes)
    element_chemical = [last(max_vp)[iid]]
    mass_chemical = [nmmi(isotopes, last(max_vp)) + msfix]
    abundance_chemical = [first(max_vp)]
    rec_addminusisotopes!(element_chemical, mass_chemical, abundance_chemical, last(max_vp), isotopes, nisotopes, iid, 2, 1, 1, first(mass_chemical), first(abundance_chemical), proportioon_cutoff, (true, true), precise)
    id = sortperm(abundance_chemical; rev = true)
    isotopes = isotopes[iid]
    (; Element = [ElementsVector(isotopes, element_chemical[i]) for i in id], Mass = mass_chemical[id], Abundance = abundance_chemical[id]) 
end

function rec_vec_ab!(els, abv, mass, tbls, abundance_cutoff, maxab, id, msn)
    msn == lastindex(tbls) && return rec_vec_ab_end!(els, abv, mass, tbls, abundance_cutoff, maxab, id)
    maxab /= first(tbls[msn].Abundance)
    pass = false
    @inbounds for (i, a) in enumerate(tbls[msn].Abundance)
        a <= 0 && break
        id[msn] = i
        next_pass = rec_vec_ab!(els, abv, mass, tbls, abundance_cutoff, maxab * a, id, msn + 1)
        pass = pass || next_pass
        next_pass ? continue : break
    end
    id[msn] = 1
    pass
end

function rec_vec_ab_end!(els, abv, mass, tbls, abundance_cutoff, maxab, id)
    maxab /= first(tbls[end].Abundance)
    pass = false
    @inbounds for (i, a) in enumerate(tbls[end].Abundance)
        ab = maxab * a
        ab < abundance_cutoff && break
        id[end] = i
        # el = [tbl.Element[i] for (tbl, i) in zip(tbls, id)]
        el = [tbl.Element[i] for (tbl, i) in zip(tbls, id)]
        push!(els, el)
        push!(abv, ab)
        push!(mass, reverse(cumsum(tbl.Mass[i] for (tbl, i) in Iterators.reverse(zip(tbls, id)))))
        pass = true
    end
    id[end] = 1
    pass
end