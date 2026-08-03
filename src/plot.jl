"""
    plot_spectrum([mz_range = nothing,] spectrum::Spectrum; deconvolution = false, abundance = 1, abtype = :max, threshold = rcrit(1e-4), kwargs...)
    plot_spectrum([mz_range = nothing,] mztable::Table; abundance = 1, abtype = :max, threshold = rcrit(1e-4), kwargs...)

Plot a spectrum.

* `mz_range::Union{Nothing, Tuple}`: nothing (indicating entire mz range) or a tuple of m/z lower bound an d upper bound.
* `deconvolution::Bool`: plot deconvoluted or convoluted spectrum.
* `abundance` sets the abundance of the peak specified by `abtype`. 
* `abtype`
    * `:max`: the largest peak.
    * `:list`: sum of listed peaks.
    * `:raw`: no abundance normalization.
* `threshold` can be a number or criteria (absolute and/or relative to maximum), data lower than this value will be set to zero. 

Other keyword arguments can controls the settings of plot.
"""
plot_spectrum(mz_range, spectrum::Spectrum; deconvolution = false, abundance = 1, abtype = :max, threshold = rcrit(1e-4), kwargs...) = 
    _plot_spectrum(mz_range, spectrum; abundance, abtype, deconvolution, threshold, kwargs...)
plot_spectrum(mz_range, mztable::Table; abundance = 1, abtype = :max, threshold = rcrit(1e-4), kwargs...) = 
    _plot_spectrum(mz_range, mztable; abundance, abtype, threshold, kwargs...)
plot_spectrum(x; kwargs...) = plot_spectrum(nothing, x; kwargs...)

"""
    plot_spectrum!([mz_range = nothing,] spectrum::Spectrum; deconvolution = false, abundance = 1, abtype = :max, threshold = rcrit(1e-4), kwargs...)
    plot_spectrum!([mz_range = nothing,] mztable::Table; abundance = 1, abtype = :max, threshold = rcrit(1e-4), kwargs...)

Plot a spectrum to an existing figure.

* `mz_range::Union{Nothing, Tuple}`: nothing (indicating entire mz range) or a tuple of m/z lower bound an d upper bound.
* `deconvolution::Bool`: plot deconvoluted or convoluted spectrum.
* `abundance` sets the abundance of the peak specified by `abtype`. 
* `abtype`
    * `:max`: the largest peak.
    * `:list`: sum of listed peaks.
    * `:raw`: no abundance normalization.
* `threshold` can be a number or criteria (absolute and/or relative to maximum), data lower than this value will be set to zero. 

Other keyword arguments can controls the settings of plot.
"""
plot_spectrum!(mz_range, spectrum::Spectrum; deconvolution = false, abundance = 1, abtype = Max(), threshold = rcrit(1e-4), kwargs...) = 
    _plot_spectrum(mz_range, spectrum; deconvolution, abundance, abtype, threshold, kwargs..., fn = plot!)
plot_spectrum!(mz_range, mztable::Table; abundance = 1, abtype = Max(), threshold = rcrit(1e-4), kwargs...) = 
    _plot_spectrum(mz_range, mztable; abundance, abtype, threshold, kwargs..., fn = plot!)
plot_spectrum!(x; kwargs...) = plot_spectrum!(nothing, x; kwargs...)

function _plot_spectrum(mz_range, mztable::Table; fn = plot, abundance = 1, abtype = Max(), threshold = rcrit(1e-4), kwargs...)
    ms, ab = plot_spectrum_params(mz_range, mztable; abundance, abtype, threshold)
    apply_plot_spectrum(ms, ab; fn, kwargs...)
end

function _plot_spectrum(mz_range, spectrum::Spectrum; fn = plot, deconvolution = false, abundance = 1, abtype = Max(), threshold = rcrit(1e-4), kwargs...)
    ms, ab = plot_spectrum_params(mz_range, spectrum; deconvolution, abundance, abtype, threshold)
    apply_plot_spectrum(ms, ab; fn, kwargs...)
end

function apply_plot_spectrum(ms, ab; fn = plot, kwargs...)
    kwargs = Dict(kwargs...)
    spec_kwargs!(kwargs)
    fn(ms, ab; kwargs...)
end
function plot_spectrum_params(mz_range, spectrum::Spectrum; deconvolution = false, abundance = 1, abtype = Max(), threshold = rcrit(1e-4), kwargs...)
    deconvolution && return plot_spectrum_params(mz_range, spectrum.table; abundance, abtype, threshold)
    if isnothing(mz_range) 
        range = length(spectrum.spectrum) * spectrum.binsize
        mz_lower = max(spectrum.initial_mass - 0.1 * range, 0)
        mz_upper = spectrum.initial_mass + 1.1 * range
    else
        mz_lower, mz_upper = mz_range 
    end
    fid, lid, delta, mzid = spectrum_id(mz_lower, mz_upper, spectrum, first(Plots.default(:size)))
    ab = map(mzid) do i 
        if i + delta < fid || i + 1 > lid
            0 
        elseif i + 1 < fid
            maximum(@view spectrum.spectrum[fid:i + delta])
        elseif i + delta > lid
            maximum(@view spectrum.spectrum[i + 1:lid])
        else
            maximum(@view spectrum.spectrum[i + 1:i + delta])
        end
    end
    abtype = abtyped(abtype)
    normalize_abundance!(ab, abundance, abtype, [Max(), List(), Raw()])
    ab_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    if isnothing(mz_range) 
        fid = findfirst(>(ab_cutoff), ab)
        lid = findlast(>(ab_cutoff), ab)
        range = ceil(Int, (lid - fid) * 0.1)
        fid = max(firstindex(ab), fid - range)
        lid = min(lastindex(ab), lid + range)
        mzid = mzid[fid:lid]
        ab = ab[fid:lid]
    end
    @. spectrum.initial_mass + (mzid * spectrum.binsize), ab
end

function plot_spectrum_params(mz_range, mztable::Table; fn = plot, abundance = 1, abtype = Max(), threshold = rcrit(1e-4))
    sp = string.(propertynames(mztable))
    colmz = lastcolnum(sp, "MZ")
    mz = getproperty(mztable, colmz)
    colab = lastcolnum(sp, "Abundance")
    if isnothing(mz_range) 
        min_mz, max_mz = extrema(mz)
        range = (max_mz - min_mz) * 0.1
        mz_lower = max(min_mz - range, 0)
        mz_upper = max_mz + range
    else
        mz_lower, mz_upper = mz_range 
    end
    id = findall(x -> mz_lower <= x <= mz_upper, mz)
    gmztable = group(getproperty(colmz), mztable[id])
    ab = map(x -> sum(getproperty(x, colab)), gmztable)
    abtype = abtyped(abtype)
    normalize_abundance!(ab, abundance, abtype, [Max(), List(), Raw()])
    ab_cutoff = minimum(makecrit_value(crit(threshold), maximum(ab)))
    mzid = findall(>=(ab_cutoff), ab)
    ab = [ab[i] for i in mzid]
    if isnothing(mz_range) 
        maxmz = maximum(mzid)
        minmz = minimum(mzid)
        range = (maxmz - minmz) / 10
        mz_lower = max(mz_lower, minmz - range)
        mz_upper = min(mz_upper, maxmz + range)
    end
    vcat(mz_lower, repeat(collect(mzid); inner = 3), mz_upper), vcat(0, [[0, a, 0] for a in ab]..., 0)
end

function spec_kwargs!(kwargs)
    get!(kwargs, :xlabel, "m/z")
    get!(kwargs, :ylabel, "Abundance")
    get!(kwargs, :title, "Spectrum")
    get!(kwargs, :label, nothing)
    kwargs
end

function spectrum_id(mz_lower, mz_upper, spectrum, np)
    mzid_s = Int((mz_lower - spectrum.initial_mass) ÷ (spectrum.binsize * spectrum.stepsize)) * spectrum.stepsize
    mzid_e = ceil(Int, (mz_upper - spectrum.initial_mass) / (spectrum.binsize * spectrum.stepsize)) * spectrum.stepsize
    fid = max(firstindex(spectrum.spectrum), mzid_s + 1)
    lid = min(lastindex(spectrum.spectrum), mzid_e + 1)
    delta, r = divrem(mzid_e - mzid_s, np)
    fid >= lid && throw(ArgumentError("No spectrum data; please adjust m/z range."))
    delta = max(delta, 1)
    fid, lid, delta, delta > 1 ? (mzid_s:delta:(mzid_e - r + delta)) : mzid_s:mzid_e
end

"""
    plot_resolving_power([mz_range = nothing,] ms::AbstractMSAnalyzer; n = 1000, label = msanalyzer_name(ms), kwargs...)

Plot the function of m/z to resolving_power.

* `mz_range::Union{Nothing, Tuple}`: nothing (indicating using `ms.mz`) or a tuple of m/z lower bound an d upper bound.
* `n::Int`: number of data points.

Other keyword arguments can controls the settings of plot.
"""
plot_resolving_power(mz_range, ms::AbstractMSAnalyzer; n = 1000, label = msanalyzer_name(ms), kwargs...) = 
    _plot_resolving_power(mz_range, ms; n, label, kwargs...)
plot_resolving_power(ms::AbstractMSAnalyzer; n = 1000, label = msanalyzer_name(ms), kwargs...) = plot_resolving_power(nothing, ms; n, label, kwargs...)

"""
    plot_resolving_power!([mz_range = nothing,] ms::AbstractMSAnalyzer; n = 1000, label = msanalyzer_name(ms), kwargs...)

Plot the function of m/z to resolving_power to an existing figure.

* `mz_range::Union{Nothing, Tuple}`: nothing (indicating using `ms.mz`) or a tuple of m/z lower bound an d upper bound.
* `n::Int`: number of data points.

Other keyword arguments can controls the settings of plot.
"""
plot_resolving_power!(mz_range, ms::AbstractMSAnalyzer; n = 1000, label = msanalyzer_name(ms), kwargs...) = 
    _plot_resolving_power(mz_range, ms; n, label, kwargs..., fn = plot!)
plot_resolving_power!(ms::AbstractMSAnalyzer; n = 1000, label = msanalyzer_name(ms), kwargs...) = plot_resolving_power!(nothing, ms; n, label, kwargs...)

function _plot_resolving_power(mz_range, ms::AbstractMSAnalyzer; fn = plot, n = 1000, label = msanalyzer_name(ms), kwargs...)
    if isnothing(mz_range)
        mz_range = ms.mz
    end
    if mz_range isa Tuple{<: Real, <: Real}
        mz_lower, mz_upper = mz_range 
    else
        throw(ArgumentError("`mz_range` must be a tuple, i.e. (lowerbound, upperbound)."))
    end
    x = mz_lower:((mz_upper - mz_lower) / n):mz_upper
    y = resolving_power.(Ref(ms), x)
    kwargs = Dict(kwargs...)
    respow_kwargs!(kwargs)
    fn(x, y; label, kwargs...)
end

function respow_kwargs!(kwargs)
    get!(kwargs, :xlabel, "m/z")
    get!(kwargs, :ylabel, "Resolution (m/Δm)")
    get!(kwargs, :title, "Resolving Power")
    kwargs
end

"""
    plot_window(window::AbstractWindow; fwhm = 0.7, binsize = 0.01, nbin_multiplier = 1, height = 0.01, label = window_name(window), kwargs...)

Plot the window function.

* `fwhm::Real`: m/z fwhm.
* `binsize::Real`: binsize of m/z value. 
* `nbin_multiplier`: multiplier of actual plotted data bin relative to `binsize`.
* `height::Real`: the lower limit of transmission. 

Other keyword arguments can controls the settings of plot.
"""
plot_window(window::AbstractWindow; fwhm = 0.7, binsize = 0.01, stepsize = 1, height = 0.01, label = window_name(window), kwargs...) = 
    _plot_window(window; fwhm, binsize, stepsize, height, label, kwargs...)

"""
    plot_window!(window::AbstractWindow; fwhm = 0.7, binsize = 0.01, nbin_multiplier = 1, height = 0.01, label = window_name(window), kwargs...)

Plot the window function to an existing figure.

* `fwhm::Real`: m/z fwhm.
* `binsize::Real`: binsize of m/z value. 
* `nbin_multiplier`: multiplier of actual plotted data bin relative to `binsize`.
* `height::Real`: the lower limit of transmission. 

Other keyword arguments can controls the settings of plot.
"""
plot_window!(window::AbstractWindow; fwhm = 0.7, binsize = 0.01, stepsize = 1, height = 0.01, label = window_name(window), kwargs...) = 
    _plot_window(window; fn = plot!, fwhm, binsize, stepsize, height, label, kwargs...)

function _plot_window(window::AbstractWindow; fn = plot, fwhm = 0.7, binsize = 0.01, stepsize = 1, height = 0.01, label = window_name(window), kwargs...)
    y = discrete_window(window, fwhm, binsize, stepsize, height)
    x = (eachindex(y) .- (length(y) ÷ 2 + 1)) .* binsize
    kwargs = Dict(kwargs...)
    window_kwargs!(kwargs)
    fn(x, y; label, kwargs...)
end

function window_kwargs!(kwargs)
    get!(kwargs, :xlabel, "m/z")
    get!(kwargs, :ylabel, "Transmission")
    get!(kwargs, :title, "Window")
    kwargs
end