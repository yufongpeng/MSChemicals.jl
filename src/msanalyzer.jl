"""
    resolving_power(ms::AbstractMSAnalyzer, mz::Real) -> Real

Resolving power of `ms` at `mz`.
"""
resolving_power(ms::AbstractMSAnalyzer, mz) = mz / fwhm_mz(ms, mz)

"""
    fwhm_mz(ms::AbstractMSAnalyzer, mz::Real) -> Real

FWHM of `mz`.
"""
fwhm_mz(ms::MSAnalyzer, mz) = ms.fwhm_mz(mz)
fwhm_mz(ms::Quadrupole, mz) = ms.fwhm 
fwhm_mz(ms::QuadrupoleIonTrap, mz) = ms.fwhm 
fwhm_mz(ms::LinearIonTrap, mz) = ms.fwhm 
fwhm_mz(ms::TOF, mz) = sqrt(mz * (3 * ms.mz50 + mz)) / ms.resolution
fwhm_mz(ms::Orbitrap, mz) = mz / ms.resolution * sqrt((3 * mz + ms.mz50 - 800) / (ms.mz50 - 200))
fwhm_mz(ms::FTICR, mz) = mz / ms.resolution * (mz + ms.mz50 - 400) / (ms.mz50 - 200)

"""
    msanalyzer_name(ms::AbstractMSAnalyzer) -> String

Common name of the MS analyzer.
"""
msanalyzer_name(::Quadrupole) = "Quadrupole"
msanalyzer_name(::QuadrupoleIonTrap) = "Quadrupole Ion Trap"
msanalyzer_name(::LinearIonTrap) = "Linear Ion Trap"
msanalyzer_name(::TOF) = "TOF"
msanalyzer_name(::Orbitrap) = "Orbitrap"
msanalyzer_name(::FTICR) = "FTICR"
msanalyzer_name(x::MSAnalyzer{T}) where T = string("MS-Analyzer with ", window_name(x.window))

measure_error(::AbstractMSAnalyzer) = [value_error, ppm_error]

msanalyzer_type(::MSAnalyzer) = MSAnalyzer
msanalyzer_type(::Quadrupole) = Quadrupole
msanalyzer_type(::QIT) = QIT
msanalyzer_type(::LIT) = LIT
msanalyzer_type(::TOF) = TOF
msanalyzer_type(::Orbitrap) = Orbitrap
msanalyzer_type(::FTICR) = FTICR

replace_mz(msanalyzer::T, mz) where {T<:AbstractMSAnalyzer} = 
    msanalyzer_type(msanalyzer)((f == :mz ? mz : getfield(msanalyzer, f) for f in fieldnames(T))...)

"""
    window_name(window::AbstractWindow) -> String

Common name of `window`.
"""
window_name(::GaussianWindow) = "Gaussian Window"
window_name(::GaussianTailedUniformWindow) = "Uniform Window (Gaussian-Tailed)"
window_name(::FixedTaperTukeyWindow) = "Tukey Window (Fixed Taper)"
window_name(::TukeyWindow) = "Tukey Window"
window_name(::CosineWindow) = "Cosine Window"
window_name(::RectWindow) = "Rectangular Window"
window_name(::PowerCosineWindow) = "Power Cosine Window"

"""
    window_parameter(msanalyzer::AbstractMSAnalyzer, mz)

Parameters of the window function `msanalyzer.window` with `mz` value. It typically returns the center and dispersion of the window. See documentation of each window for details.
"""
window_parameter(msanalyzer::AbstractMSAnalyzer{<:GaussianTailedUniformWindow}, μ) = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{<:GaussianWindow}, μ) = (μ, fwhm_mz(msanalyzer, μ) ^ 2 / log(256))
window_parameter(msanalyzer::AbstractMSAnalyzer{<:FixedTaperTukeyWindow}, μ) = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{<:TukeyWindow}, μ) = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{<:CosineWindow}, μ) = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{<:PowerCosineWindow}, μ) = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{<:RectWindow}, μ) = (μ, fwhm_mz(msanalyzer, μ))
window_parameter(msanalyzer::AbstractMSAnalyzer{<:SuperGaussianWindow}, μ) = (μ, fwhm_mz(msanalyzer, μ) / (sqrt(8) * log(2) ^ (1 / msanalyzer.window.power)))

(::SampledGaussianWindow)(ν, t) = sgaussian_norm(ν, t)
(::GaussianWindow)(x, μ, σ2) = gaussian_norm(x, μ, σ2)
(f::SuperGaussianWindow)(x, μ, σ) = supergaussian_norm(x, μ, σ, f.power)
function (f::FixedTaperTukeyWindow)(x, μ::S, fwhm::T) where {S, T}
    taperwidth = min(f.taperwidth, fwhm)
    p = abs(x - μ)
    plateau = fwhm - taperwidth
    if p > plateau / 2 + taperwidth
        zero(promote_type(float(S), float(T)))
    elseif p > plateau / 2
        (cos(p * π / taperwidth) + 1) / 2
    else
        one(promote_type(float(S), float(T)))
    end
end

function (f::TukeyWindow)(x, μ::S, fwhm::T) where {S, T}
    width = fwhm / (2 - f.taperproportion)
    plateau = width * (1 - f.taperproportion)
    p = abs(x - μ) 
    if p > width
        zero(promote_type(float(S), float(T)))
    elseif p > plateau
        taperwidth = width - plateau
        p = p - plateau
        (cos(p * π / taperwidth) + 1) / 2
    else
        one(promote_type(float(S), float(T)))
    end
end

(::CosineWindow)(x, μ::S, fwhm::T) where {S, T} = abs(x - μ) > fwhm ? zero(promote_type(float(S), float(T))) : (cos(π * abs(x - μ) / fwhm) + 1) / 2
(::RectWindow)(x, μ, fwhm) = abs(x - μ) < fwhm / 2 ? 1 : 0
function (f::PowerCosineWindow)(x, μ::S, fwhm::T) where {S, T}
    p = (x - μ) * 2 / fwhm * acos(2 ^ (1 - 1 / f.power) - 1)
    # p = (x - μ) * 2 / fwhm
    abs(p) > π ? zero(promote_type(float(S), float(T))) : (0.5 + cos(p) / 2) ^ f.power
end
    
function (f::GaussianTailedUniformWindow)(x, μ::S, fwhm::T) where {S, T}
    p = abs(x - μ) 
    p == 0 && return one(promote_type(float(S), float(T)))
    taperwidth = min(f.taperwidth, fwhm)
    plateau = fwhm - taperwidth
    if plateau == 0
        gaussian_norm(p, 0, fwhm ^ 2 / log(256))
    elseif p > plateau / 2
        gaussian_norm(p - plateau / 2, 0, taperwidth ^ 2 / log(256))
    else
        one(promote_type(float(S), float(T)))
    end
end

gaussian_norm(x::Real, μ::Real, σ2::Real) = exp(-(x - μ) ^ 2 / 2 / σ2)
supergaussian_norm(x::Real, μ::Real, σ::Real, n::Real) = exp(-((x - μ) / sqrt(2) / σ) ^ n)
sgaussian_norm(ν::Real, t) = exp(-ν ^ 2 / 2 / t)
sgaussian_norm(ν::AbstractRange, t::T) where T = sgaussian_norm!(zeros(T, length(ν)), ν, t)
function sgaussian_norm!(u::AbstractVector, ν::AbstractRange, t) 
    delta = step(ν)
    start = ν.start
    num = exp(-start ^ 2 / 2t)
    deltanum1 = exp(-delta ^ 2 / 2t)
    deltanum2 = exp(-delta / t)
    prev_ν = start - delta
    prev_num = num / deltanum1 / deltanum2 ^ prev_ν
    for i in eachindex(ν)
        u[i] = prev_num * deltanum1 * deltanum2 ^ prev_ν
        prev_num = u[i]
        prev_ν += delta
    end
    u
end

"""
    discrete_window(msanalyzer::AbstractMSAnalyzer, mz::Vector, binsize::Real, nbin_multiplier::Real, height::Real) 
    discrete_window(window::AbstractWindow, fwhm::Real, binsize::Real, nbin_multiplier::Real, height::Real)

Discrete values of `window` with `fwhm` or ms analyzer `msanalyzer.window` at `mz` values (FWHM is computed through function [`fwhm_mz`](@ref)). 
`binsize` sets the size of m/z bin, `nbin_multiplier` sets the interval between sampled bins, and `height` sets minimal window value.
"""
function discrete_window(msanalyzer::AbstractMSAnalyzer, μ::Vector{T}, binsize, nbin_multiplier, height) where T
    fwhm_ref = fwhm_mz(msanalyzer, first(μ))
    window_ref = discrete_window(msanalyzer.window, fwhm_ref, binsize, nbin_multiplier, height)
    windows = Vector{Vector{T}}(undef, length(μ))
    for (i, m) in enumerate(μ)
        fwhm = fwhm_mz(msanalyzer, m)
        if fwhm - fwhm_ref < binsize * nbin_multiplier
            windows[i] = window_ref 
        else
            windows[i] = discrete_window(msanalyzer.window, fwhm, binsize, nbin_multiplier, height)
            window_ref = windows[i]
            fwhm_ref = fwhm
        end
    end
    windows
end

function discrete_window(window::CosineWindow, fwhm, binsize, nbin_multiplier, height)
    dfwtm = ceil(fwhm * 2 / binsize)
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
    k = [window(x * binsize, 0, fwhm) for x in 0:dfwtm]
    vcat(reverse(k), k[begin + 1:end])
end

function discrete_window(window::PowerCosineWindow, fwhm, binsize, nbin_multiplier, height)
    dfwtm = ceil(fwhm * 2 / binsize)
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
    k = [window(x * binsize, 0, fwhm) for x in 0:dfwtm]
    vcat(reverse(k), k[begin + 1:end])
end

function discrete_window(window::TukeyWindow, fwhm, binsize, nbin_multiplier, height)
    dfwtm = ceil(fwhm / (2 - window.taperproportion) / binsize)
    plateau = ceil(dfwtm * (1 - window.taperproportion))
    cosfwhm = dfwtm - plateau
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
    k = [CosineWindow()(x, 0, cosfwhm) for x in 0:dfwtm - plateau]
    vcat(reverse(k), [1 for _ in 1:(2 * plateau)], k[begin + 1:end])
end

function discrete_window(window::RectWindow, fwhm, binsize, nbin_multiplier, height)
    width = ceil(fwhm / 2 / binsize)
    p, r = divrem(width, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : width
    k = [Int(x <= width) for x in 0:dfwtm]
    vcat(reverse(k), k[begin + 1:end])
end

function discrete_window(window::FixedTaperTukeyWindow, fwhm, binsize, nbin_multiplier, height)
    dfwtm = ceil((fwhm + window.taperwidth) / 2 / binsize)
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
    plateau = fwhm - window.taperwidth
    k = [window(x * binsize, 0, plateau) for x in 0:dfwtm]
    vcat(reverse(k), k[begin + 1:end])
end

function discrete_window(window::GaussianTailedUniformWindow, fwhm, binsize, nbin_multiplier, height)
    taperstart = ceil(Int, (fwhm - window.taperwidth) / 2 / binsize)
    offset = rem(taperstart, nbin_multiplier)
    fwhm = fwhm - taperstart * binsize * 2
    k = _discrete_window(SampledGaussianWindow(), fwhm, binsize, nbin_multiplier, height; offset)
    vcat(reverse(k), ones(2 * taperstart), k[begin + 1:end])
end

discrete_window(window::GaussianWindow, fwhm, binsize, nbin_multiplier, height) = discrete_window(SampledGaussianWindow(), fwhm, binsize, nbin_multiplier, height)
function discrete_window(window::SuperGaussianWindow, fwhm, binsize, nbin_multiplier, height)
    hwhm = fwhm / 2
    bsqrtt = hwhm / binsize
    ihwhml = floor(Int, hwhm / binsize)
    ihwhmr = ceil(Int, hwhm / binsize)
    ihwhmr = ihwhmr > ihwhml ? ihwhmr : ihwhml + 1
    t = estimate_t(window, bsqrtt, ihwhml, ihwhmr)
    dfwtm = round(Int, sqrt(2) * t * log(1 / height) ^ (1 / window.power)) + 1
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r > 0 ? (p + 1) * nbin_multiplier : dfwtm
    k = [window(x, 0, t) for x in 0:dfwtm]
    vcat(reverse(k), k[begin + 1:end])
end
function _discrete_window(window::SampledGaussianWindow, fwhm, binsize, nbin_multiplier, height; offset = 0)
    hwhm = fwhm / 2
    bsqrtt = hwhm / binsize
    ihwhml = floor(Int, hwhm / binsize)
    ihwhmr = ceil(Int, hwhm / binsize)
    ihwhmr = ihwhmr > ihwhml ? ihwhmr : ihwhml + 1
    t = estimate_t(window, bsqrtt, ihwhml, ihwhmr)
    dfwtm = round(Int, sqrt(2t * log(1 / height))) + 1
    p, r = divrem(dfwtm, nbin_multiplier)
    dfwtm = r == offset ? dfwtm : r > offset ? (p + 1) * nbin_multiplier + offset : p * nbin_multiplier + offset
    window(0:dfwtm, t)
end

function discrete_window(window::SampledGaussianWindow, fwhm, binsize, nbin_multiplier, height)
    k = _discrete_window(window, fwhm, binsize, nbin_multiplier, height)
    vcat(reverse(k), k[begin + 1:end])
end

function estimate_t(window::SampledGaussianWindow, bsqrtt, ihwhml, ihwhmr, log_n = 0)
    t = bsqrtt ^ 2 / log(4)
    log_n > 1000 && return t
    # p = window(0, t) / 2
    hl = window(ihwhml, t)
    hr = window(ihwhmr, t)
    if hr <= 0.5 <= hl 
        return t
    elseif 0.5 > hl
        bsqrtt += (0.5 - (hl + hr) / 2) / (hl - hr)
        return estimate_t(window, bsqrtt, ihwhml, ihwhmr, log_n + 1)
    else
        bsqrtt += (0.5 - (hl + hr) / 2) / (hl - hr)
        return estimate_t(window, bsqrtt, ihwhml, ihwhmr, log_n + 1)
    end
end

function estimate_t(window::SuperGaussianWindow, bsqrtt, ihwhml, ihwhmr, log_n = 0)
    t = bsqrtt / sqrt(2) / log(2) ^ (1 / window.power)
    log_n > 1000 && return t
    # p = window(0, 0, t) / 2
    hl = window(ihwhml, 0, t)
    hr = window(ihwhmr, 0, t)
    if hr <= 0.5 <= hl 
        return t
    elseif 0.5 > hl
        bsqrtt += (0.5 - (hl + hr) / 2) / (hl - hr)
        return estimate_t(window, bsqrtt, ihwhml, ihwhmr, log_n + 1)
    else
        bsqrtt += (0.5 - (hl + hr) / 2) / (hl - hr)
        return estimate_t(window, bsqrtt, ihwhml, ihwhmr, log_n + 1)
    end
end