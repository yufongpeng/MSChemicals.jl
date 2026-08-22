"""
    coelutingisobars_filter! = ms_filter! ∘ elution_filter!

"""
coelutingisobars_filter!(ci::CoelutingIsobars) = ms_filter!(elution_filter!(ci))

"""
    elution_filter!(ci::CoelutingIsobars)

Filter `ci.isobars` for each target chemicals using `ci.elution` and store the result in `ci.tables`.
"""
function elution_filter!(ci::CoelutingIsobars)
    exp = ci.target
    lib = ci.isobar
    exptransitions = seriesanalyzedchemical.(ci.target.Chemical)
    libtransitions = seriesanalyzedchemical.(ci.isobar.Chemical)
    # [elution, target/isobar]
    exp_elu_value = [[fn(first(c)) for c in exptransitions] for (fn, _) in ci.elution]
    lib_elu_value = [[fn(first(c)) for c in libtransitions] for (fn, _) in ci.elution]
    # [target, isobar]
    quals = _elution_filter!(ci.elution, exp_elu_value, lib_elu_value)
    tables = empty!(ci.tables)
    for (i, qual) in enumerate(quals)
        j = findall(qual)
        filter!(k -> !ischemicalequal(ci.target.Chemical[i], ci.isobar.Chemical[k]), j)
        push!(tables, ci.isobar[j])
    end
    ci
end

function _elution_filter!(
        elution, 
        exp_elu_value, 
        lib_elu_value, 
        quals = [trues(length(first(lib_elu_value))) for i in eachindex(first(exp_elu_value))]
    )
    all(x -> all(!, x), quals) && return quals
    i = 0
    while i < lastindex(elution)
        i += 1 
        i1 = min(lastindex(exp_elu_value), i)
        i2 = min(lastindex(lib_elu_value), i)
        for j in eachindex(quals)
            all(!, quals[j]) && continue 
            _ielution_filter!(
                last(elution[i]), 
                exp_elu_value[i1][j],
                lib_elu_value[i2],
                quals[j]
            )
        end
        any(any, quals) && continue 
    end
    quals
end

function _ielution_filter!(crit_elu::Criteria, target_elu, lib_elu, quals)
    tol = union(makecrit_delta(crit_elu, target_elu)...)
    quals .*= [x in tol for x in lib_elu]
end

"""
    ms_filter!(ci::CoelutingIsobars)

Filter `ci.tables` for each target chemicals using `ci.msanalyzer`. Filtering only occurs on transitions that corresponding msanalyzer has a vector of target mz (this supposes to be empty, but it's fine with some elements).  
"""
function ms_filter!(ci::CoelutingIsobars)
    (isempty(ci.tables) || length(ci.target) != length(ci.tables)) && throw(ArgumentError("Run `elution_filter!` first."))
    msanalyzer = first.(ci.msanalyzer)
    ab_cutoff = last.(ci.msanalyzer)
    ti = findall(x -> x.mz isa Vector, msanalyzer)
    cmz = [mz.(seriesanalyzedchemical(x)) for x in ci.target.Chemical]
    for i in eachindex(ci.target)
        _msanalyzer = deepcopy(msanalyzer)
        mztable = ci.tables[i]
        for j in ti 
            push!(_msanalyzer[j].mz, cmz[i][j])
        end
        for (j, ma) in enumerate(_msanalyzer)
            if j in ti 
                mztable = Isolation(ma, mztable; stage = j, threshold = ab_cutoff[j])
            else
                mztable = AllIons(mztable)
            end
        end
        ci.tables[i] = mztable 
    end
    ci
end

"""
    isobar_table(ci::CoelutingIsobars; isotope = "[13C]", threshold = rcrit(1e-4), error_elution = measure_error.(ci.elution), error_mz = measure_error.(ci.msanalyzer)) -> Table

A table listing all target chemicals and their isobars grouped by isomeric state based on `isotope` from `ci`.

# Arguments
* `isotope::AbstractString`: a minor isotope.
* `threshold` can be a number or criteria, representing the lower limit of proportion of isobar affecting the target chemical. 
* `error_elution::Vector{<:Vector}`: error function(s) of each elution dimension.
* `error_mz::Vector{<:Vector}`: error function(s) of each msanalyzer.
"""
function isobar_table(ci::CoelutingIsobars; 
        isotope = "[13C]",
        threshold = rcrit(1e-4), 
        error_elution = measure_error.(first.(ci.elution)), 
        error_mz = measure_error.(first.(ci.msanalyzer))
    )
    id = findall(!isempty, ci.tables)
    # group by parent and In 
    gt = [group_isotopologues(mztable; isotope) for mztable in ci.tables[id]]
    cols = Pair{Symbol, Vector}[]
    push!(cols, :Chemical => vcat((repeat([ci.target.Chemical[i]], length(gt[j])) for (j, i) in enumerate(id))...))
    target_mz = vcat((repeat([mz.(seriesanalyzedchemical(ci.target.Chemical[i]))], length(gt[j])) for (j, i) in enumerate(id))...)
    push!(cols, :Isobar => vcat(getproperty.(gt, :Chemical)...))
    for (i, fns) in enumerate(error_elution)
        elu_value = [(first(ci.elution[i])(x), first(ci.elution[i])(y)) for (x, y) in zip(last(first(cols)), last(last(cols)))]
        for fn in fns
            push!(cols, Symbol(measure_name(first(ci.elution[i]), fn)) => [fn(x...) for x in elu_value])
        end
    end
    sp = string.(propertynames(first(gt)))
    colmz = allcolnum(sp, "MZ")
    colab = allcolnum(sp, "Abundance")
    for (i, fns) in enumerate(error_mz)
        mz_value = zip([x[i] for x in target_mz], vcat(getproperty.(gt, colmz[i])...))
        for fn in fns
            push!(cols, Symbol(measure_name(colmz[i], fn)) => [fn(x...) for x in mz_value])
        end
    end
    for c in colab
        push!(cols, Symbol(string(c, "(%)")) => vcat(getproperty.(gt, c) ./ getproperty(ci.target[id], c) .* 100...))
    end
    table = Table(; cols...)
    isnothing(threshold) && return table
    abundance_cutoff = minimum(makecrit_value(crit(threshold), 1)) * 100
    table[findall(>=(abundance_cutoff), getproperty(table, first(last(cols))))]
end

# lex(v1, v2, i = 1) = i > lastindex(v1) || (i <= lastindex(v2) && (v1[i] == v2[i] ? lex(v1, v2, i + 1) : v1[i] < v2[i])) 
# elution dim: indivisual fwhm / crit
# [rt => Criteria]
# [rt => Resolution(fwhm) => cutoff]
# ms: analyzer dependent / crit
# [msanalyzer => Criteria (abundance)], [Criteria]
# error, percentage error, ppm error, resolution