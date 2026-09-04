"""
    Ionization(chemicaltable::Table; [adduct, abundance, proportion,] adduction = AdductIon, chemicalparser = ChemicalExpressionParser(), adductparser = AdductParser(), threshold = rcrit(1e-4), threading = nothing, kwargs...) -> Table
    Ionization(chemical::AbstractString; chemicalparser = ChemicalExpressionParser(), kwargs...) -> Table
    Ionization(chemical::AbstractChemical; adduct, adduction = AdductIon, adductparser = AdductParser(), kwargs...) -> Table
    Ionization(chemical::AbstractAdductIon; kwargs...) -> Table

Ionization of `chemicaltable.Chemical`. 
The returned table contains isotopologues of adduct ions and can be further processed by other spectrum related function such as [`MSScan`](@ref), [`Isolation`](@ref), and etc. 

# Keyword Arguments
* `adduction`: constructor for adduct ion. 
* `chemicalparser::AbstractChemicalParser`: parser for `string` chemical input. The default parser is `ChemicalTransitionParser(ChemicalExpressionParser())` (See [`ChemicalTransitionParser`](@ref)).
* `adductparser::AbstractAdductParser`: parser for `string` adduct input. The default parser is `AdductParser()` (See [`AdductParser`](@ref)).
* `adduct`: adducts being parsed by `adductparser`. The parsed adduct is a named tuple and used as keyword arguments of function `ionize` with argument `adduction` and each chemical. For individualizing adducts, use column `Adduct` in `chemicaltable`.
* `abundance` sets the abundance of the chemical. It can also be column `Abundance` in `chemicaltable`. 
* `proportion`: proportion of adduct ion relative to original chemical. For individualizing adducts, use column `Proportion` in `chemicaltable`. The length of each elements should macthes to that of `adduct`.
* `threshold` can be a number or [`Criteria`](@ref), representing the lower limit of abundance (absolute and/or relative to maximal value of each spectrum). 
* `threading`: force to use multiple threads (`true`) or single thread (`false`); `nothing` lets the program determine. 
Other `kwargs` are used by function [`Isotopologues`](@ref).
"""
function Ionization(mztable::Table; adduction = AdductIon, threading = nothing, chemicalparser = ChemicalExpressionParser(), adductparser = AdductParser(), threshold = rcrit(1e-4), kwargs...)
    :Chemical in propertynames(mztable) || throw(ArgumentError("No column `Chemical` in input table."))
    chemical = eltype(mztable.Chemical) <: AbstractChemical ? mztable.Chemical : parse_chemical.(Ref(chemicalparser), mztable.Chemical; kwargs...)
    kwargs = Dict(kwargs...)
    if :Adduct in propertynames(mztable)
        adduct = map(vectorize, mztable.Adduct)
        delete!(kwargs, :adduct)
    elseif haskey(kwargs, :adduct)
        adduct = vectorize(kwargs[:adduct], length(mztable))
        delete!(kwargs, :adduct)
    else
        throw(ArgumentError("Require adducts information. Use column `Adduct` of the table or keyword arguments `adduct`."))
    end
    adduct = map(x -> parse_adduct.(Ref(adductparser), vectorize(x)), adduct)
    if :Abundance in propertynames(mztable)
        abundance = mztable.Abundance
        delete!(kwargs, :abundance)
    elseif haskey(kwargs, :abundance)
        abundance = vectorize(kwargs[:abundance], length(mztable))
        delete!(kwargs, :abundance)
    else
        abundance = [1.0 for _ in eachindex(mztable)]
    end
    if :Proportion in propertynames(mztable)
        proportion = mztable.Proportion
        delete!(kwargs, :proportion)
    elseif haskey(kwargs, :proportion)
        proportion = vectorize(kwargs[:proportion], length(mztable))
        delete!(kwargs, :proportion)
    else
        proportion = [repeat([1.0], length(first(adduct))) for _ in eachindex(mztable)]
    end
    proportion = map(vectorize, proportion)
    id = vcat(([(i, j) for j in eachindex(adduct[i])] for i in eachindex(mztable))...)
    rn = min(length(id), Threads.nthreads())
    if rn == 1 
        threading = false 
    end
    if isnothing(threading)
        ab = mean(abundance)
        s = mean(length(chemicalelements(x)) for x in mztable.Chemical)
        b = mean(mean(last, chemicalelements(x)) for x in mztable.Chemical)
        b = sum(b ^ (1/x) for x in 1:msstage(first(mztable.Chemical)))
        # println((0.02b + 0.2sqrt(-2 * b * log2(min(1, minimum(makecrit_value(crit(threshold), ab)) / ab)))) ^ 1.2s * (rn - 1))
        threading = (0.02b + 0.2sqrt(-2b * log2(min(1, minimum(makecrit_value(crit(threshold), ab)) / ab)))) ^ 1.2s * (rn - 1) > 1e6
    end
    if threading
        t = Vector{Table}(undef, length(id))
        Threads.@threads for k in eachindex(t)
            i, j = id[k]
            t[k] = Isotopologues(ionize(adduction, chemical[i]; adduct[i][j]...); kwargs..., id = (k, ), abundance = abundance[i] * proportion[i][j], threshold, iter = true, sort = true)
        end
        tbl = Table(; (p => ChainedVector(getproperty.(t, p)) for p in propertynames(t[1]))...)
    else
        t = [Isotopologues(ionize(adduction, chemical[i]; adduct[i][j]...); kwargs..., id = (k, ), abundance = abundance[i] * proportion[i][j], threshold, iter = true, sort = true) for (k, (i, j)) in enumerate(id)]
        tbl = Table(; (p => ChainedVector(getproperty.(t, p)) for p in propertynames(t[1]))...)
    end
    colab = lastcolnum(propertynames(tbl), "Abundance"; error = false)
    ab = getproperty(tbl, colab)
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>=(abundance_cutoff), ab)
    tbl[id]
end

Ionization(v::Vector; kwargs...) = Ionization(Table(; Chemical = v); kwargs...) 
Ionization(x::AbstractString; chemicalparser = ChemicalExpressionParser(), kwargs...) = Ionization(parse_chemical(chemicalparser, x); kwargs...) 
Ionization(x::AbstractChemical; adduction = AdductIon, adduct, adductparser = AdductParser(), kwargs...) = Isotopologues(ionize(adduction, x; parse_adduct(adductparser, adduct)...); kwargs...) 
Ionization(x::AbstractAdductIon; kwargs...) = Isotopologues(x; kwargs...) 
Ionization(::Isobars; kwargs...) = throw(ArgumentError("`Isobars` is not supported by `Ionization`."))
Ionization(::Isotopomers; kwargs...) = throw(ArgumentError("`Isotopomers` is not supported by `Ionization`."))

"""
    MSScan(msanalyzer = TOF(), mztable; min_bin_fwhm = 50) -> Union{Spectrum, Nothing}
    MSScan(msanalyzer = TOF(), spectrum; min_bin_fwhm = 50) -> Union{Spectrum, Nothing}

Perform MS scan. Theoretical m/z signals are convoluted with `msanalyzer.window`. 

# Arguments
* `msanalyzer::AbstractMSAnalyzer`: a MS analyzer to perform MS Scan. See documentation of specific analyzer for detailed settings.
* `mztable::Table`: a table containing columns
    * `Abundance1`, `Abundance2`, ..., `Abundancen`. The last column will be utilized.
    * `MZ1`, `MZ2`, ..., `MZn`. The last column will be utilized.
* `spectrum::Spectrum`: `spectrum.table` is utilized as `mztable`.
* `min_bin_fwhm::Int`: minimal number of bins within fwhm. 
"""
function MSScan(msanalyzer::AbstractMSAnalyzer, mztable::Table; min_bin_fwhm = 50)
    isempty(mztable) && return nothing
    digits = hasproperty(msanalyzer, :digits) ? msanalyzer.digits : nothing
    stepsize = hasproperty(msanalyzer, :stepsize) ? msanalyzer.stepsize : nothing
    if !isnothing(stepsize) && !isnothing(digits) 
        unit = round(10.0 ^ (-digits); digits)
        p = round(Int, stepsize / unit)
        isapprox(stepsize - p * unit, 0) || throw(ArgumentError("`stepsize` must be mutiples of `10 ^ (-digits)`."))
    end
    sp = string.(propertynames(mztable))
    colmz = lastcolnum(sp, "MZ")
    colab = lastcolnum(sp, "Abundance")
    id = sortperm(getproperty(mztable, colmz))
    if !isnothing(msanalyzer.mz)
        lower_mz, upper_mz = msanalyzer.mz
        mz_vector = getproperty(mztable, colmz)
        filter!(x -> lower_mz < mz_vector[x] < upper_mz, id)
    end
    mztable = mztable[id]
    mz_vector = getproperty(mztable, colmz)
    ab_vector = getproperty(mztable, colab)
    fwhm = fwhm_mz(msanalyzer, first(mz_vector))
    binsize = min(makecrit_value(crit(msanalyzer.accuracy), first(mz_vector))..., fwhm / min_bin_fwhm)
    nbin_multiplier = 1
    if !isnothing(stepsize)
        nbin_multiplier = ceil(Int, stepsize / binsize)
        binsize = stepsize / nbin_multiplier
    end
    kernels = discrete_window(msanalyzer, mz_vector, binsize, nbin_multiplier, minimum(ab_vector) / maximum(ab_vector) / 10)
    # ihwhm, k = convolution_window(window, hwhm, binsize, nbin_multiplier, minimum(mztable.Abundance) / 10; normalize = true)
    initmass = isnothing(digits) ? first(mz_vector) : round(first(mz_vector); digits)
    binmass, ibins = binnify(mz_vector, binsize, initmass)
    spectrum_offset = bin_offset(initmass, first(binmass), binsize, nbin_multiplier)
    convolution_vector = [a .* k for (a, k) in zip(ab_vector, kernels)]
    ibins .+= Int(length(first(kernels))) ÷ 2 + 1
    spectrum = zeros(eltype(ab_vector), last(ibins) + Int(length(last(kernels))) ÷ 2)
    @inbounds for (c, i) in zip(convolution_vector, ibins)
        convolution_offset = Int(length(c) ÷ 2)
        spectrum[i - convolution_offset:i + convolution_offset] .+= c 
    end
    # spectrum = spectrum_offset > 1 ? spectrum[spectrum_offset:end] : spectrum
    spectrum = spectrum_offset > 0 ? vcat(zeros(spectrum_offset), spectrum) : spectrum
    # ibins .-= spectrum_offset - 1
    ibins .+= spectrum_offset
    mztable = Table(mztable; 
        Bin_MZ = binmass, 
        Convolution = convolution_vector, 
        Bin_id = ibins
    )
    initial_mass = first(binmass) - (first(ibins) - 1) * binsize
    # pt = peak_table(mztable, spectrum, initial_mass, binsize)
    Spectrum(spectrum, initial_mass, binsize, nbin_multiplier, mztable)
end
MSScan(mztable; min_bin_fwhm = 50) = MSScan(TOF(), mztable; min_bin_fwhm)
MSScan(msanalyzer::AbstractMSAnalyzer, spec::Spectrum; min_bin_fwhm = 50) = MSScan(msanalyzer, spec.table; min_bin_fwhm)

"""
    AllIons(mz_range = nothing, mztable) -> Table
    AllIons(mz_range = nothing, spectrum) -> Table

Allow all Ions within m/z range entering the next MS stage. 

# Arguments
* `mz_range::Union{Nothing, Tuple}`: nothing (indicating all ions) or a tuple of m/z lower bound and upper bound.
* `mztable::Table`: a table containing columns
    * `MZ1`, `MZ2`, ..., `MZn`. The last column will be utilized.
* `spectrum::Spectrum`: `spectrum.table` is utilized as `mztable`.
"""
function AllIons(mz_range, mztable::Table)
    colmz = lastcolnum(propertynames(mztable), "MZ")
    if !isnothing(mz_range)
        lower_mz, upper_mz = mz_range 
        id = findall(x -> lower_mz < x < upper_mz, getproperty(mztable, colmz))
        mztable = mztable[id]
    end
    mztable
end
AllIons(mztable) = AllIons(nothing, mztable)
AllIons(mz_range, spec::Spectrum) = AllIons(mz_range, spec.table)

"""
    Isolation(msanalyzer, mztable; stage = nothing, threshold = rcrit(1e-4)) -> Table
    Isolation(msanalyzer, spectrum; stage = nothing, threshold = rcrit(1e-4)) -> Table

Isolating target ion(s) with specific m/z values and resolutions to enter the next MS stage. 

# Arguments
* `msanalyzer::AbstractMSAnalyzer`: a MS analyzer to perform MS filtering. See documentation of specific analyzer for detailed settings.
* `mztable::Table`: a table containing columns    
    * `ID`: ID tuples. Each elements represents ID number of ions of each MS stage. 
    * `Abundance1`, `Abundance2`, ..., `Abundancen`. 
    * `MZ1`, `MZ2`, ..., `MZn`. 
* `spectrum::Spectrum`: `spectrum.table` is utilized as `mztable`.
* `stage::Int`: MS stage. Default `nothing` for the last MS stage.
* `threshold` can be a number or [`Criteria`](@ref) (absolute and/or relative to maximum), representing the lower limit of abundance. 
"""
function Isolation(msanalyzer::AbstractMSAnalyzer, mztable::Table; stage = nothing, threshold = rcrit(1e-4)) 
    isempty(mztable) && return mztable
    sp = string.(propertynames(mztable))
    if isnothing(stage)
        colmz = lastcolnum(sp, "MZ")
        colab = lastcolnum(sp, "Abundance")
    else
        colmz = ithcolnum(sp, "MZ", stage)
        colab = ithcolnum(sp, "Abundance", stage)
    end
    _Isolation(msanalyzer, mztable, colmz, colab, threshold)
end
Isolation(msanalyzer::AbstractMSAnalyzer, spec::Spectrum; stage = nothing, threshold = rcrit(1e-4)) = Isolation(msanalyzer, spec.table; stage, threshold)

function _Isolation(msanalyzer::AbstractMSAnalyzer, mztable::Table, colmz::Symbol, colab::Symbol, threshold)
    params = [window_parameter(msanalyzer, mz) for mz in vectorize(msanalyzer.mz)]
    isempty(params) && return mztable
    isempty(mztable) && return mztable
    ab = map(mztable) do r
        getproperty(r, colab) * maximum(msanalyzer.window(getproperty(r, colmz), param...) for param in params)
    end
    tab = minimum(makecrit_value(crit(threshold), maximum(ab)))
    id = findall(>(tab), ab)
    Table(mztable; [colab => ab]...)[id]
end

"""
    SelectedIonMonitor(transitiontable, mztable; threading = nothing, threshold = rcrit(1e-4)) -> Table
    SelectedIonMonitor(transitiontable, spectrum; threading = nothing, threshold = rcrit(1e-4)) -> Table

Selected ion monitoring. 

# Arguments
* `transitiontable::Table`: each row represents a transition. Use column `Transition` (optional) for specifying transition name. 
Other columns must be in analysis-fragmentation-analysis order. 
Analysis columns contain MS analyzers and Fragmentation columns contain producttables (See `Fragmentation` for detail).
* `mztable::Table`: a table containing columns    
    * `ID`: ID tuples. Each elements represents ID number of ions of each MS stage. 
    * `Abundance1`
    * `MZ1`
* `spectrum::Spectrum`: `spectrum.table` is utilized as `mztable`.
* `threshold` can be a number or [`Criteria`](@ref) (absolute and/or relative to maximum), representing the lower limit of abundance. 
* `threading`: force to use multiple threads (`true`) or single thread (`false`); `nothing` lets the program determine. 
"""
function SelectedIonMonitor(transitiontable::Table, mztable::Table; threading = nothing, threshold = rcrit(1e-4)) 
    isempty(mztable) && return transitiontable
    sp = string.(propertynames(mztable))
    colmz = lastcolnum(sp, "MZ")
    colmz == :MZ1 || throw(ArgumentError("Require only single MZ column `MZ1."))
    colab = lastcolnum(sp, "Abundance")
    colab == :Abundance1 || throw(ArgumentError("Require only single Abundance column `Abundance1."))
    transitions = Table(transitiontable; Transition = nothing)
    tables = _SelectedIonMonitor(columns(transitions), mztable, colmz, colab, threshold, threading)
    Table(transitiontable; MZTable = tables)
end
SelectedIonMonitor(transitiontable::Table, spec::Spectrum; threshold = rcrit(1e-4)) = SelectedIonMonitor(transitiontable, spec.table; threshold)

function _SelectedIonMonitor(transitions, mztable::Table, colmz::Symbol, colab::Symbol, threshold, threading)
    length(transitions) % 2 < 1 && throw(ArgumentError("Transitions must have 'analysis-fragmentation-analysis...` pattern."))
    isolation = true 
    mztable = repeat([mztable], length(first(transitions)))
    for transition in transitions
        if isolation
            mztable = [_Isolation(trans, mzt, colmz, colab, threshold) for (trans, mzt) in zip(transition, mztable)]
        else
            mztable = [Fragmentation(trans, mzt; threading, threshold) for (trans, mzt) in zip(transition, mztable)]
            colmz = Symbol(string("MZ", parse(Int, string(match(r"\d+", string(colmz)).match)) + 1))
            colab = Symbol(string("Abundance", parse(Int, string(match(r"\d+", string(colab)).match)) + 1))
        end
        isolation = !isolation
    end
    mztable
end

"""
    Fragmentation(product_table, precursor_table; threading = nothing, threshold = rcrit(1e-4)) -> Table
    Fragmentation(product_table, spectrum; threading = nothing, threshold = rcrit(1e-4)) -> Table

Fragmentation of `precursor_table.Chemical` or `spectrum.table.Chemical` into `product_table.Product`.

# Arguments
* `producttable::Table`: a table containing columns
    * `Product`: products of each precursor. 
    * `Proportion`: proportion of fragmentation of each product relative to precursor signal. This column is optional; the default is that each product share precursor signals equally. 
* `precursortable::Table`: a table containing columns
    * `ID`: ID tuples. Each element represents ID number of ions of each MS stage. 
    * `Chemical`: chemical objects. 
    * `Abundance1`, `Abundance2`, ..., `Abundancen`. The last column will be utilized.
    * `MZ1`, `MZ2`, ..., `MZn`. The last column will be utilized.
* `spectrum::Spectrum`: `spectrum.table` is utilized as `precursortable`.
* `threshold` can be a number or [`Criteria`](@ref), representing the lower limit of abundance (absolute and/or relative to maximal value of each spectrum). 
* `threading`: force to use multiple threads (`true`) or single thread (`false`); `nothing` lets the program determine. 
"""
function Fragmentation(producttable::Table, mztable::Table; chemicalparser = ChemicalExpressionParser(), threading = nothing, threshold = rcrit(1e-4))
    # group -> precursor_table, elements_precursor
    :ID in propertynames(mztable) || throw(ArgumentError("No column `ID` in precursor_table."))
    :Chemical in propertynames(mztable) || throw(ArgumentError("No column `Chemical` in precursor_table."))
    :Abundance1 in propertynames(mztable) || throw(ArgumentError("No column `Abundance1`, ..., `Abundacnen` in precursor_table."))
    :MZ1 in propertynames(mztable) || throw(ArgumentError("No column `MZ1`, ..., `MZn` in precursor_table."))
    if isempty(mztable) 
        lc = lastcolnum(propertynames(mztable), "MZ")
        (n, ) = match(r"MZ(\d+)", string(lc))
        n = parse(Int, n)
        mspre = map(1:n) do i 
            s = Symbol(string("MZ", i))
            s => getproperty(mztable, s)
        end
        abpre = map(1:n) do i 
            s = Symbol(string("Abundance", i))
            s => getproperty(mztable, s)
        end
        return Table(; ID = mztable.ID, Chemical = mztable.Chemical, mspre..., [Symbol(string("MZ", n + 1)) => getproperty(mztable, :MZ1)]..., abpre..., [Symbol(string("Abundance", n + 1)) => getproperty(mztable, :Abundance1)]...) 
    end
    in(:Chemical, propertynames(producttable)) || throw(ArgumentError("No column `Chemical` in product_table."))
    if !(eltype(producttable.Chemical) <: AbstractChemical)
        producttable = Table(producttable; Chemical = [parse_chemical(chemicalparser, x) for x in producttable.Chemical])
    end
    :Product in propertynames(producttable) || throw(ArgumentError("No column `Product` in product_table."))
    if any(x -> !(eltype(x) <: AbstractChemical), producttable.Product)
        producttable = Table(producttable; Product = [[parse_chemical(chemicalparser, y) for y in x] for x in producttable.Product])
    end
    allequal(msstage, mztable.Chemical) || throw(ArgumentError("Chemicals have to be in the same MS stage."))
    all(x -> all(y -> msstage(y) < 2, x), producttable.Product) || throw(ArgumentError("Products should not be MS/MS pairs."))
    if !in(:Proportion, propertynames(producttable))
        producttable = Table(producttable; Proportion = [1 for _ in eachindex(producttable)])
    end
    # MS1 threshold
    # threshold = acrit(minimum(makecrit_value(crit(threshold), maximum(mztable.Abundance1))))
    gmztable = group(getproperty(:ID), mztable)
    prodid = map(gmztable) do mzt
        findfirst(x -> ischemicalequal(x, detectedchemical(chemicalparent(first(mzt.Chemical)))), producttable.Chemical)
    end
    ks = findall(!isnothing, prodid)
    rn = min(length(ks), Threads.nthreads())
    if isnothing(threading)
        colab = lastcolnum(propertynames(mztable), "Abundance")
        ab = mean(getproperty(mztable, colab))
        s = mean(length(chemicalelements(first(gmztable[x].Chemical))) for x in ks)
        b = mean(mean(last, chemicalelements(first(gmztable[x].Chemical))) for x in ks)
        b = sum(b ^ (1/x) for x in 1:(msstage(first(mztable.Chemical)) + 1))
        np = mean(length(producttable.Product[prodid[x]]) for x in ks)
        # println((0.02b + 0.2sqrt(-2 * b * log2(min(1, minimum(makecrit_value(crit(threshold), ab)) / ab)))) ^ 1.2s * (rn - 1))
        threading = np * (0.02b + 0.2sqrt(-2b * log2(min(1, minimum(makecrit_value(crit(threshold), ab)) / ab)))) ^ 1.2s * (rn - 1) > 1e6
    end
    if threading
        t = Vector{Table}(undef, length(ks))
        ks = collect(ks)
        Threads.@threads for i in eachindex(t)
            precursor_table = gmztable[ks[i]]
            pid = prodid[ks[i]]
            t[i] = TandemIsotopologues(chemicalparent(first(precursor_table.Chemical)); 
                threshold,
                precursor_table,
                product = producttable.Product[pid],
                proportion = producttable.Proportion[pid],
                transmission = sum(producttable.Proportion[pid]),
                iter = true,
                sort = true
            )
        end
        Table(; (p => ChainedVector(getproperty.(t, p)) for p in propertynames(first(t)))...)
    else
        t = Vector{Table}(undef, length(ks))
        for (i, k) in enumerate(ks)
            precursor_table = gmztable[k]
            pid = prodid[k]
            t[i] = TandemIsotopologues(chemicalparent(first(precursor_table.Chemical)); 
                threshold,
                precursor_table,
                product = producttable.Product[pid],
                proportion = producttable.Proportion[pid],
                transmission = sum(producttable.Proportion[pid]),
                iter = true,
                sort = true
            )
        end
        Table(; (p => ChainedVector(getproperty.(t, p)) for p in propertynames(first(t)))...)
    end
end
Fragmentation(producttable::Table, spec::Spectrum; kwargs...) = Fragmentation(producttable, spec.table; kwargs...)

"""
    peak_table(spectrum::Spectrum; alg = LocalMaxima(), abundance = 1, abtype = :max, threshold = rcrit(1e-4)) -> Table
    peak_table(transitiontable::Table; groupedisotopomers = true, isotope = "[13C]") -> Table

Extract peaks from a [`Spectrum`](@ref) or [`SelectedIonMonitor`](@ref).

# Arguments
* `abundance` sets the abundance of the peak specified by `abtype`. 
* `abtype`.
    * `:max`: the largest peak.
    * `:list`: sum of listed peaks.
    * `:raw`: no abundance normalization.
* `threshold` can be a number or [`Criteria`](@ref) (absolute and/or relative to `abundance`), representing the lower limit of abundance. 
* `groupedisotopomers`: whether group isotopologues by isotopomer state based on `isotope`.
* `isotope::String`: minor isotope.
"""
peak_table(spectrum::Spectrum; alg = LocalMaxima(), abundance = 1, abtype = Max(), threshold = rcrit(1e-4)) = peak_table(spectrum.table, spectrum.spectrum, spectrum.initial_mass, spectrum.binsize, spectrum.stepsize; alg, abundance, abtype, threshold)
function peak_table(mztable::Table, spectrum, initial_mass, binsize, nbin_multiplier; alg = LocalMaxima(), abundance = 1, abtype = Max(), threshold = rcrit(1e-4))
    maxbin = [find_nearest_peak(alg, spectrum, r.Bin_id, r.Convolution) for r in mztable] 
    if nbin_multiplier > 1
        maxbin = map(maxbin) do mbin 
            isnothing(mbin) && return nothing
            r = rem(mbin - 1, nbin_multiplier)
            nbin_multiplier - r < r ? mbin + nbin_multiplier - r : mbin - r 
        end
    end
    id = findall(!isnothing, maxbin)
    mztable = Table(mztable[id]; Max_bin = maxbin[id])
    c = length(first(mztable.Convolution)) ÷ 2 + 1
    gmztable = group(getproperty(:Max_bin), mztable)
    tuples = map(pairs(gmztable)) do (ibin, smztable) 
        cab = map(smztable) do r 
            r.Convolution[c + ibin - r.Bin_id]
        end
        (; Chemical = Isobars(getproperty(smztable, :Chemical), cab), MZ = initial_mass + (ibin - 1) * binsize, Abundance = spectrum[ibin])
    end
    table = Table([t for t in tuples])
    abtype = abtyped(abtype)
    normalize_abundance!(table.Abundance, abundance, abtype, [Max(), List(), Raw()])
    abundance_cutoff = minimum(makecrit_value(crit(threshold), maximum(table.Abundance)))
    table[findall(>=(abundance_cutoff), table.Abundance)]
end

function peak_table(transitiontable::Table; groupedisotopomers = true, isotope = "[13C]")
    :MZTable in propertynames(transitiontable) || throw(ArgumentError("No column `MZTable` in transitiontable"))
    id = findall(!isempty, transitiontable.MZTable)
    mztables = @view transitiontable.MZTable[id]
    if groupedisotopomers
        tables = map(mztables) do table 
            group_isotopologues(table; isotope)
        end
    else
        tables = mztables
    end
    sp = string.(propertynames(first(tables)))
    colmz = allcolnum(sp, "MZ")
    colab = allcolnum(sp, "Abundance")
    gt = map(tables) do mztable 
        isobar = Isobars(mztable.Chemical, hcat(getproperty.(Ref(mztable), colab)...))
        transitions = chemicaltransition.(mztable.Chemical)
        if groupedisotopomers
            (; Chemical = isobar, [c => mean(getproperty(mztable, c), weights(getproperty(mztable, d))) for (c, d) in zip(colmz, colab)]..., [c => sum(getproperty(mztable, c)) for c in colab]...)
        else
            gid, gt = unique_group_mz_ab(mztable, transitions, colmz, colab)
            (; Chemical = isobar, first(gt)...)
        end
    end
    if :Transition in propertynames(transitiontable)
        Table(Table(; Transition = transitiontable.Transition[id]), Table(collect(NamedTuple, gt)))
    else
        Table(collect(NamedTuple, gt))
    end
end

function find_nearest_peak(alg::LocalMaxima, convolution, i, k)
    j = findfirst(>(alg.threshold * maximum(k)), k)
    ihwhm = floor(Int, length(k) / 2) - j + 1
    peak = [convolution[i], convolution[i]]
    dir = [true, true]
    start = [false, false]
    ibin = [i, i]
    for j in 0:ihwhm
        if first(dir)
            if convolution[i - j - 1] > first(peak)
                start[begin] = true
                peak[begin] = convolution[i - j - 1]
            else
                dir[begin] = false
                ibin[begin] = i - j
            end
        end
        if last(dir)
            if convolution[i + j + 1] > last(peak) 
                start[end] = true
                peak[end] = convolution[i + j + 1]
            else
                dir[end] = false
                ibin[end] = i + j
            end
        end
        dir'start > 0 || break 
    end
    r = @. (!)(dir) * start
    id = if !any(start)
        ibin[begin]
    elseif all(r) && peak[begin] < peak[end]
        ibin[end]
    elseif all(r)
        ibin[begin]
    elseif first(r)
        ibin[begin]
    elseif last(r)
        ibin[end]
    else
        nothing 
    end
    id
end

function bin_offset(outmass, binmass, binsize, nbin_multiplier)
    x = round(Int, (outmass - binmass) / binsize)
    # x == 0 ? 1 : x > 0 ? x + 1 : nbin_multiplier + x + 1
    x == 0 ? x : x > 0 ? nbin_multiplier - x : -x
end

function binnify(mass, binsize, init = first(mass))
    binmass = similar(mass)
    ibins = zeros(Int, length(binmass))
    ibin = 0
    for (i, m) in enumerate(mass)
        while init < m 
            ibin += 1
            init += binsize
        end
        if init - m > binsize / 2 
            binmass[i] = init - binsize 
            ibins[i] = ibin - 1
        else
            binmass[i] = init
            ibins[i] = ibin
        end
    end
    binmass, ibins .- first(ibins)
end