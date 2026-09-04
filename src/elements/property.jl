"""
    iselement(x::AbstractString) -> Bool

Determine if `x` is an element.
"""
iselement(x::AbstractString) = haskey(elements_isotopes(), x)

"""
    isisotope(x::AbstractString) -> Bool

Determine if `x` is an isotope (including element).
"""
isisotope(x::AbstractString) = haskey(elements_mass(), x)

"""
    ismajor(x::AbstractString) -> Bool

Determine if `x` is the major isotope.
"""
ismajor(x::AbstractString) = x == major_isotope(x)

"""
    isminor(x::AbstractString, i::Int = 1) -> Bool

Determine if `x` is the `i`th minor isotope.
"""
isminor(x::AbstractString, i::Int = 1) = x == minor_isotope(x, i)

"""
    parent_element(x::AbstractString) -> String

Parent elements of isotope `x`.
"""
parent_element(x::AbstractString) = get(elements_parents(), x, "")

"""
    major_isotope(x::AbstractString) -> String

Major isotope of isotope `x`.
"""
function major_isotope(x::AbstractString) 
    e = parent_element(x)
    if haskey(elements_isotopes(), e)
        first(elements_isotopes()[e])
    else
        "" 
    end
end

"""
    minor_isotope(x::AbstractString, i = 1) -> String

`i`th minor isotope of isotope `x`.
"""
function minor_isotope(x::AbstractString, i::Int = 1)
    e = parent_element(x)
    if haskey(elements_isotopes(), e)
        v = elements_isotopes()[e]
        i < lastindex(v) ? v[i + 1] : "" 
    else
        ""
    end
end