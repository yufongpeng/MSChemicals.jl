"""
    Criteria{A, R}

A type for creating criteria of numeric value. Criterion may be a cutoff value, allowed tolerance relative to an input value, or acceptable range (interval(s)).

# Fields
* `aval`: criterion for the absolute value (unnormalized).
* `rval`: criterion for the relative value (normalized by an input)

# Examples
```julia
julia> peak_crit = Criteria(10, 0.2);

julia> qualified_peak(x, basepeak, crit) = all(<=(x), makecrit_value(crit, basepeak)); # filter out small peaks based on base peak. 

julia> makecrit_value(peak_crit, 100) # create criteria for peak value 100
(10, 20)

julia> makecrit_value(peak_crit, 10) # create criteria for peak value 10
(10, 2)

julia> qualified_peak(10, 100, peak_crit)
false

julia> qualified_peak(20, 100, peak_crit)
true

julia> qualified_peak(10, 10, peak_crit)
true

julia> qualified_peak(2, 10, peak_crit)
false

julia> qc_peak(x, refpeak, crit) = all(in(x), makecrit_delta(crit, refpeak)); # filter out deviated peaks.

julia> makecrit_delta(peak_crit, 100) # create criteria for reference peak value 100
(Intervals.Interval{Float64, Intervals.Closed, Intervals.Closed}(90.0, 110.0), Intervals.Interval{Float64, Intervals.Closed, Intervals.Closed}(80.0, 120.0))

julia> qc_peak(121, 100, peak_crit)
false

"""
struct Criteria{A, R}
    aval::A
    rval::R
end

abstract type AbstractAbundance end 
struct Input <: AbstractAbundance end
struct Max <: AbstractAbundance end
struct List <: AbstractAbundance end 
struct Total <: AbstractAbundance end 
struct Raw <: AbstractAbundance end 