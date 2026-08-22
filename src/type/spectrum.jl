"""
    Spectrum{T}

Spectrum type.

# Fields
* `spectrum::Vector{T}`: spectrum intensity vector.
* `initial_mass::T`: the first m/z value of `spectrum`.
* `binsize::T`: the size of m/z bin, i.e. m/z difference of consecutive elements.
* `stepsize::Int`: number of bins for each data point of the actual spectrum. 
* `table::Table`: the source table of ions. 
"""
struct Spectrum{T}
    spectrum::Vector{T}
    initial_mass::T 
    binsize::T
    stepsize::Int
    table::Table
end

# Record procedure
# struct MSAnalysis 
#     procedure::Vector
#     table::Table
#     tables::Vector
# end

abstract type AbstractWindow end
abstract type ZeroCenteredtWindow end

"""
    GaussianWindow <: AbstractWindow 

Gaussian window.

# Window Function
    (::GaussianWindow)(x, μ, σ2)

## Arguments
* `μ`: mean (center).
* `σ2`: variance.
"""
struct GaussianWindow <: AbstractWindow end

struct GaussianTailedUniformWindow{T} <: AbstractWindow
    taperwidth::T
end

"""
    FixedTaperTukeyWindow <: AbstractWindow 

Tukey window with fixed width of taper.

# Fields 
* `taperwidth`: width of taper.

# Window Function
    (::FixedTaperTukeyWindow)(x, μ, fwhm)

## Arguments
* `μ`: mean (center).
* `fwhm`: full width at half maximum.
"""
struct FixedTaperTukeyWindow{T} <: AbstractWindow 
    taperwidth::T
end

"""
    TukeyWindow <: AbstractWindow 

Tukey window.

# Window Function
    (::TukeyWindow)(x, μ, fwhm)

## Arguments
* `μ`: mean (center).
* `fwhm`: full width at half maximum.
"""
struct TukeyWindow{T} <: AbstractWindow 
    taperproportion::T
end

"""
    CosineWindow <: AbstractWindow 

Cosine window.

# Window Function
    (::CosineWindow)(x, μ, fwhm)

## Arguments
* `μ`: mean (center).
* `fwhm`: full width at half maximum.
"""
struct CosineWindow <: AbstractWindow end

"""
    RectWindow <: AbstractWindow 

Rectangular window.

# Window Function
    (::RectWindow)(x, μ, fwhm)

## Arguments
* `μ`: mean (center).
* `fwhm`: full width at half maximum.
"""
struct RectWindow <: AbstractWindow end

"""
    SampledGaussianWindow <: ZeroCenteredtWindow 

Sampled gaussian window.

# Window Function
    (::SampledGaussianWindow)(x, t)

## Arguments
* `t`: variance.
"""
struct SampledGaussianWindow <: ZeroCenteredtWindow end

"""
    SuperGaussianWindow <: AbstractWindow 

Super gaussian window.

# Fields 
* `power`: the power of exponential. `2` is equivalent to typical gaussian window.

# Window Function
    (::SuperGaussianWindow)(x, μ, σ)

## Arguments
* `μ`: mean (center).
* `σ`: standard deviation.
"""
struct SuperGaussianWindow{T} <: AbstractWindow 
    power::T
end

"""
    PowerCosineWindow{T} <: AbstractWindow 

Power cosine (top half) window.

# Fields 
* `power`: the power of cosine.

# Window Function
    (::PowerCosineWindow)(x, μ, fwhm)

## Arguments
* `μ`: mean (center).
* `fwhm`: full width at half maximum.
"""
struct PowerCosineWindow{T} <: AbstractWindow 
    power::T
end

abstract type AbstractMSAnalyzer{W, M} end 

"""
    MSAnalyzer{W, M} <: AbstractMSAnalyzer{W, M} 

Generic MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s). It can be a 
    * `Real`: center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: m/z range.
    * `Nothing`: m/z value not specified.
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `fwhm_mz`: a function mapping m/z value to fwhm.
* `digits`: digits of m/z value. `nothing` indicates no limit on digits.
* `stepsize`: step size for discrete mass scan such as quadrupole. It must be multiples of `10 ^ (-digits)`. If `digits` is `nothing`, `stepsize` is `nothing` as well.
"""
struct MSAnalyzer{W, M} <: AbstractMSAnalyzer{W, M}
    window::W
    mz::M
    accuracy
    fwhm_mz 
    digits
    stepsize
end

"""
    Quadrupole{W, M} <: AbstractMSAnalyzer{W, M} 

Quadrupole MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s). It can be a 
    * `Real`: center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: m/z range.
    * `Nothing`: m/z value not specified.
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `fwhm::Real`: m/z fwhm.
* `digits`: digits of m/z value. `nothing` indicates no limit on digits.
* `stepsize`: step size for mass scan. It must be multiples of `10 ^ (-digits)`. If `digits` is `nothing`, `stepsize` is `nothing` as well.

# Constructors
    Quadrupole(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, digits = nothing, stepsize = nothing, unit = 0.7, flatness = 100, taperproportion = 0.2)

## Arguments
* `offset::Real`: m/z offset of the center of window. The resulting m/z center will be `mz + offset`.
* `unit::Real`: unit mass resolution (fwhm); the largest fwhm that remains gaussian. Window function becomes super gaussian (flat top) for fwhm larger than unit.
* `flatness::Real`: how flat is the gaussian fall-off of the super gaussian function. Practical values range from 10 to 1000. 
* `taperproportion::Real`: the proportion of the window that is tapered for `TukeyWindow` (used in single m/z isolation). Practical values range from 0.1 to 8. 
"""
struct Quadrupole{W, M} <: AbstractMSAnalyzer{W, M}
    window::W
    mz::M
    accuracy
    fwhm
    digits
    stepsize
end 

function Quadrupole(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, digits = nothing, stepsize = nothing, unit = 0.7, flatness = 100, taperproportion = 0.2)
    mz = isnothing(mz) ? mz : @.(mz + offset)
    if isnothing(mz) || mz isa Tuple
        if fwhm > unit
            f = sqrt(log(2, flatness))
            n = log(2, f) / log(2, unit / fwhm * (f - 1) + 1) * 2
            Quadrupole(SuperGaussianWindow(n), mz, accuracy, fwhm, digits, stepsize)
        else
            Quadrupole(GaussianWindow(), mz, accuracy, fwhm, digits, stepsize)
        end
    elseif taperproportion >= 1
        Quadrupole(CosineWindow(), mz, accuracy, fwhm, digits, stepsize)
    elseif taperproportion > 0
        Quadrupole(TukeyWindow(taperproportion), mz, accuracy, fwhm, digits, stepsize)
    else
        Quadrupole(RectWindow(), mz, accuracy, fwhm, digits, stepsize)
    end 
end

doc_qit = """
    QuadrupoleIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    const QIT = QuadrupoleIonTrap

Quadrupole Ion Trap (QIT) MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s). It can be a 
    * `Real`: center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: m/z range.
    * `Nothing`: m/z value not specified.
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `fwhm::Real`: m/z fwhm.
* `digits`: digits of m/z value. `nothing` indicates no limit on digits.
* `stepsize`: step size for mass scan. It must be multiples of `10 ^ (-digits)`. If `digits` is `nothing`, `stepsize` is `nothing` as well.

# Constructors
    QuadrupoleIonTrap(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, power = 4)

## Arguments
* `offset::Real`: m/z offset of the center of window. The resulting m/z center will be `mz + offset`.
* `power::Real`: the power of super gaussian function (used in SWIFT isolation). Practical values range from 4 to 8. 
"""
@doc doc_qit
struct QuadrupoleIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    fwhm
    digits
    stepsize
end 

function QuadrupoleIonTrap(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, power = 4)
    mz = isnothing(mz) ? mz : @.(mz + offset)
    if isnothing(mz) || mz isa Tuple
        # Instability scan mode
        QuadrupoleIonTrap(GaussianWindow(), mz, accuracy, fwhm, digits, stepsize)
    else
        # SWIFT
        QuadrupoleIonTrap(SuperGaussianWindow(power), mz, accuracy, fwhm, digits, stepsize)
    end
end
@doc doc_qit
const QIT = QuadrupoleIonTrap

doc_lit = """
    LinearIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    const LIT = LinearIonTrap

Linear Ion Trap (LIT) MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s). It can be a 
    * `Real`: center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: m/z range.
    * `Nothing`: m/z value not specified.
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `fwhm::Real`: m/z fwhm.
* `digits`: digits of m/z value. `nothing` indicates no limit on digits.
* `stepsize`: step size for mass scan. It must be multiples of `10 ^ (-digits)`. If `digits` is `nothing`, `stepsize` is `nothing` as well.

# Constructors
    LinearIonTrap(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, taperproportion = 0.2)

## Arguments
* `offset::Real`: m/z offset of the center of window. The resulting m/z center will be `mz + offset`.
* `taperproportion::Real`: the proportion of the window that is tapered for `TukeyWindow` (used in SWIFT isolation).
"""
@doc doc_lit
struct LinearIonTrap{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    fwhm
    digits
    stepsize
end 

function LinearIonTrap(mz = nothing; accuracy = 0.1, offset = 0, fwhm = 0.7, stepsize = nothing, taperproportion = 0.2)
    mz = isnothing(mz) ? mz : @.(mz + offset)
    if isnothing(mz) || mz isa Tuple
        # Instability scan mode
        LinearIonTrap(GaussianWindow(), mz, accuracy, fwhm, digits, stepsize)
    elseif taperproportion >= 1
        LinearIonTrap(CosineWindow(), mz, accuracy, fwhm, digits, stepsize)
    elseif taperproportion > 0
        LinearIonTrap(TukeyWindow(taperproportion), mz, accuracy, fwhm, digits, stepsize)
    else
        LinearIonTrap(RectWindow(), mz, accuracy, fwhm, digits, stepsize)
    end 
end
@doc doc_lit
const LIT = LinearIonTrap

doc_tof = """
    TimeOfFlight{W, M} <: AbstractMSAnalyzer{W, M} 
    const TOF = TimeofFlight

Time-of-Flight (TOF) MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s). It can be a 
    * `Real`: center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: m/z range.
    * `Nothing`: m/z value not specified.
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `resolution`: theoretical maximal m/z resolution. 
* `mz50`: m/z value with half `resolution`.

# Constructors
    TimeOfFlight(mz = nothing; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 80000, mz50 = 120)
"""
@doc doc_tof
struct TimeOfFlight{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    resolution
    mz50
end 
TimeOfFlight(mz = nothing; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 80000, mz50 = 120) = TimeOfFlight(window, mz, accuracy, resolution, mz50)
# 100000 .* sqrt.(x) ./ (8 .+ sqrt.(x)) .- 20000
@doc doc_tof
const TOF = TimeOfFlight

"""
    Orbitrap{W, M} <: AbstractMSAnalyzer{W, M} 

Orbitrap MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s). It can be a 
    * `Real`: center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: m/z range.
    * `Nothing`: m/z value not specified.
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `resolution`: m/z resolution @200. 
* `mz50`: m/z value with half `resolution`.

# Constructors
    Orbitrap(mz = nothing; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 480000, mz50 = 800)
"""
struct Orbitrap{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    resolution
    mz50
end 
Orbitrap(mz = nothing; window = GaussianWindow(), accuracy = rcrit(1e-6), resolution = 480000, mz50 = 800) = Orbitrap(window, mz, accuracy, resolution, mz50)

doc_fticr = """
    FourierTransformIonCyclotronResonance{W, M} <: AbstractMSAnalyzer{W, M} 
    const FTICR = FourierTransformIonCyclotronResonance

Fourier-transform ion cyclotron resonance (FTICR) MS Analyzer.

# Fields
* `window::AbstractWindow`: window function.
* `mz`: analyzed m/z value(s). It can be a 
    * `Real`: center m/z value.
    * `Vector{<:Real}`: multiple center m/z values.
    * `Tuple{<:Real, <:Real}`: m/z range.
    * `Nothing`: m/z value not specified.
* `accuracy`: m/z accuracy, either a number or a `Criteria`.
* `resolution`: m/z resolution @200. 
* `mz50`: m/z value with half `resolution`.

# Constructors
    FourierTransformIonCyclotronResonance(mz = nothing; window = GaussianWindow(), accuracy = rcrit(2e-7), resolution = 20000000, mz50 = 400)
"""
@doc doc_fticr
struct FourierTransformIonCyclotronResonance{W, M} <: AbstractMSAnalyzer{W, M} 
    window::W
    mz::M
    accuracy
    resolution
    mz50
end 
FourierTransformIonCyclotronResonance(mz = nothing; window = GaussianWindow(), accuracy = rcrit(2e-7), resolution = 20000000, mz50 = 400) = FourierTransformIonCyclotronResonance(window, mz, accuracy, resolution, mz50)
@doc doc_fticr
const FTICR = FourierTransformIonCyclotronResonance

abstract type AlgSpecPeak end
struct LocalMaxima <: AlgSpecPeak
    threshold::Real
end 
LocalMaxima() = LocalMaxima(0.1)

"""
    CoelutingIsobars

A type containing settings and results of coelution isobars. 

# Fields
* `elution::Vector{<:Pair}`: elution function (e.g. `retentiontime`)-criteria pairs. 
* `msanalyzer::Vector{<:Pair}`: msanalyzer-abundance criteria pairs.
* `target::Table`: a table of target chemicals in the form of `Isotopomers`.
* `isobar::Table`: a table of potential isobars in the form of `Isotopomers`.
* `tables::Vector`: tables of isotopologues of each target chemicals.

The index of `msanalyzer` matches the index of transition of chemicals. Filtering only occurs on transitions that corresponding msanalyzer has a vector of target mz (this supposes to be empty, but it's fine with some elements). If `msanalyzer` is shorter, no filtering is applied for the transitions that have no corresponding msanalyzer.

# Constructors
    CoelutingIsobars(elution, msanalyzer, table::Table; ci_filter = 2, kwargs...)
    CoelutingIsobars(elution, msanalyzer, target::Table, isobar::Table; ci_filter = 2)

## Arguments
* `table::Table`: table of target chemicals not in the form of `Isotopomers`; `target` and `isobar` are generated from this table. 
* `ci_filter::Int`: number of step to filter isobars. `0` means no filtering, and `tables` would be empty; `1` means filtering by only `elution`; `2` means filtering by both `elution` and `msanalyzer`.
* Other keyword arguments are for `Isotopologues` to generate `isobar` from `table`.
"""
struct CoelutingIsobars 
    elution::Vector{<:Pair}
    msanalyzer::Vector{<:Pair}
    target::Table 
    isobar::Table
    tables::Vector
end

function CoelutingIsobars(elution, ms, table::Table; ci_filter = 2, kwargs...) 
    isobar = TandemIsotopologues(table; kwargs...)
    target = isobar[[findfirst(x -> ischemicalequal(x, y), isobar.Chemical) for y in table.Chemical]]
    CoelutingIsobars(elution, ms, target, isobar; ci_filter)
end

function CoelutingIsobars(elution, msanalyzer, target::Table, isobar::Table; ci_filter = 2) 
    ci = CoelutingIsobars(elution, msanalyzer, target, isobar, Any[])
    ci_filter < 1 && return ci 
    elution_filter!(ci)
    ci_filter -= 1
    ci_filter < 1 && return ci 
    ms_filter!(ci)
end