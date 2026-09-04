
"""
    unique_elements(elements::Vector{<:Pair}) -> Vector{<:Pair}
    unique_elements(elements::Dict) -> Dict
    unique_elements(elements::Dictionary) -> Dictionary
    unique_elements(elements::ElementsVector) -> ElementsVector

Elements container with no duplicated element keys or zeros.
"""
unique_elements(elements::T) where T = unique_elements(T, elements)
unique_elements(::Type{<:Dictionary}, elements) = filter!(!=(0), dictionary_elements(Dictionary, elements))
unique_elements(::Type{<:Dict}, elements) = Dict(pairs(unique_elements(Dictionary, elements)))
unique_elements(::Type{<:Vector{<:Pair}}, elements) = collect(pairs(unique_elements(Dictionary, elements)))
unique_elements(::Type{<:Dict}, elements::Dict) = filter(x -> last(x) != 0, elements)
unique_elements(::Type{<:Vector{<:Pair}}, elements::Dict) = filter!(x -> last(x) != 0, collect(elements))
unique_elements(::Type{<:Dictionary}, elements::Dictionary) = filter(!=(0), elements)
function unique_elements(::Type{<:ElementsVector}, elements::ElementsVector) 
    id = findall(!=(0), elements.numbers)
    ElementsVector(elements.elements[id], elements.numbers[id])
end
unique_elements(::Type{<:Vector{<:Pair}}, elements::ElementsVector) = filter!(x -> last(x) != 0, collect(elements))

sort_unique_elements(x) = sort!(unique_elements(x))

"""
    dictionary_elements(Dicttype = Dict, elements::Vector{<:Pair}) -> Dicttype
    dictionary_elements(Dicttype = Dict, elements::Dict) -> Dicttype
    dictionary_elements(Dicttype = Dict, elements::Dictionary) -> Dicttype
    dictionary_elements(Dicttype = Dict, elements::ElementsVector) -> Dicttype

Create a dictionary from any elements containers. As elements can be duplicated in some containers and no duplication in a dictionary, the new dictionary is convenient for updating elements number.
"""
dictionary_elements(elements) = dictionary_elements(Dict, elements)
dictionary_elements(::Type{Dict}, elements::Vector{<:Pair}) = Dict(pairs(dictionary_elements(Dictionary, elements)))
dictionary_elements(::Type{Dictionary}, elements::Vector{<:Pair}) = groupsum(first, last, elements)
dictionary_elements(::Type{Dict}, elements::Dict) = elements
dictionary_elements(::Type{Dictionary}, elements::Dict) = Dictionary(keys(elements), values(elements))
dictionary_elements(::Type{Dict}, elements::Dictionary) = Dict(pairs(elements))
dictionary_elements(::Type{Dictionary}, elements::Dictionary) = elements
dictionary_elements(::Type{Dict}, elements::ElementsVector) = Dict(elements)
dictionary_elements(::Type{Dictionary}, elements::ElementsVector) = Dictionary(elements.elements, elements.numbers)

"""
    gain_elements(elements::Vector{<:Pair}, y...) -> Vector{<:Pair}
    gain_elements(elements::Dict, y...) -> Dict
    gain_elements(elements::Dictionary, y...) -> Dictionary
    gain_elements(elements::ElementVector, y...) -> ElementVector

Add elements in `y` to copied `elements`.
"""
gain_elements(elements, y...) = gain_elements!(copy(elements), y...)

"""
    gain_elements!(elements::Vector{<:Pair}, y...) -> Vector{<:Pair}
    gain_elements!(elements::Dict, y...) -> Dict
    gain_elements!(elements::Dictionary, y...) -> Dictionary
    gain_elements!(elements::ElementVector, y...) -> ElementVector

Add elements in `y` to `elements`.
"""
function gain_elements!(elements::Dict, y...) 
    for d in y 
        _gain_elements!(elements, d)
    end
    filter!(!=(0), elements)
end

function gain_elements!(elements::Dictionary, y...) 
    for d in y 
        _gain_elements!(elements, d)
    end
    filter!(!=(0), elements)
end

function gain_elements!(elements::Vector{<:Pair}, y...)
    for d in y 
        _gain_elements!(elements, d)
    end
    elements
end

function gain_elements!(elements::ElementsVector, y...)
    for d in y 
        _gain_elements!(elements, d)
    end
    elements
end

_gain_elements!(elements, y::Dictionary) = __gain_elements!(elements, pairs(y))
_gain_elements!(elements, y) = __gain_elements!(elements, y)

function __gain_elements!(elements::Dict, y)
    for (k, v) in y
        get!(elements, k, 0)
        elements[k] += v
    end
end

function __gain_elements!(elements::Dictionary, y)
    for (k, v) in y
        get!(elements, k, 0)
        elements[k] += v
    end
end

function __gain_elements!(elements::Vector{<:Pair}, y)
    for k in y
        last(k) != 0 && push!(elements, k)
    end
end

function __gain_elements!(elements::ElementsVector, y)
    for (k, v) in y
        i = findfirst(==(k), elements.elements)
        if isnothing(i)
            push!(elements.elements, k)
            push!(elements.numbers, v)
        else
            elements.numbers[i] = elements.numbers[i] + v 
        end
    end
    elements
end

function parallel_gain_elements!(els::Vector{Vector{ElementsVector}})
    prev_el = String[]
    prev_ns = Vector{Int}[] 
    @inbounds for i in Iterators.reverse(eachindex(els[begin]))
        if isempty(prev_el) 
            prev_el = els[begin][i].elements
            prev_ns = [el[i].numbers for el in els]
            continue
        end
        iter_el = els[begin][i].elements
        curr_el = els[begin][i].elements
        curr_ns = [copy(el[i].numbers) for el in els]
        pushed = false
        for (ie, e) in enumerate(prev_el) 
            j = findfirst(==(e), iter_el)
            if isnothing(j)
                if pushed
                    push!(curr_el, e)
                else
                    curr_el = push!(copy(curr_el), e)
                    pushed = true
                end
                for (curr_n, prev_n) in zip(curr_ns, prev_ns)
                    push!(curr_n, prev_n[ie])
                end
            else 
                for (curr_n, prev_n) in zip(curr_ns, prev_ns)
                    curr_n[j] = curr_n[j] + prev_n[ie]
                end
            end
        end
        for (j, el) in enumerate(els) 
            el[i] = ElementsVector(curr_el, curr_ns[j])
        end
        prev_el = curr_el
        prev_ns = curr_ns
    end
    els
end

"""
    loss_elements(elements::Vector{<:Pair}, y...) -> Vector{<:Pair}
    loss_elements(elements::Dict, y...) -> Dict
    loss_elements(elements::Dictionary, y...) -> Dictionary
    loss_elements(elements::ElementVector, y...) -> ElementVector

Substract elements in `y` from copied `elements`.
"""
loss_elements(elements, y...) = loss_elements!(copy(elements), y...)

"""
    loss_elements!(elements::Vector{<:Pair}, y...) -> Vector{<:Pair}
    loss_elements!(elements::Dict, y...) -> Dict
    loss_elements!(elements::Dictionary, y...) -> Dictionary
    loss_elements!(elements::ElementVector, y...) -> ElementVector

Substract elements in `y` from `elements`.
"""
function loss_elements!(elements::Dict, y...) 
    for d in y 
        _loss_elements!(elements, d)
    end
    filter!(!=(0), elements)
end

function loss_elements!(elements::Dictionary, y...) 
    for d in y 
        _loss_elements!(elements, d)
    end
    filter!(!=(0), elements)
end

function loss_elements!(elements::Vector{<:Pair}, y...)
    for d in y 
        _loss_elements!(elements, d)
    end
    elements
end

function loss_elements!(elements::ElementsVector, y...)
    for d in y 
        _loss_elements!(elements, d)
    end
    elements
end

_loss_elements!(elements, y::Dictionary) = __loss_elements!(elements, pairs(y))
_loss_elements!(elements, y) = __loss_elements!(elements, y)

function __loss_elements!(elements::Dict, y)
    for (k, v) in y
        get!(elements, k, 0)
        elements[k] -= v
    end
end

function __loss_elements!(elements::Dictionary, y)
    for (k, v) in y
        get!(elements, k, 0)
        elements[k] -= v
    end
end

function __loss_elements!(elements::Vector{<:Pair}, y)
    for (k, v) in y
        v != 0 && push!(elements, k => -v)
    end
end

function __loss_elements!(elements::ElementsVector, y)
    for (k, v) in y
        i = findfirst(==(k), elements.elements)
        if isnothing(i)
            push!(elements.elements, k)
            push!(elements.numbers, v)
        else
            elements.numbers[i] = elements.numbers[i] - v 
        end
    end
    elements
end

function universal_loss_elements!(els::Vector{Vector{ElementsVector}})
    prev_el = els[begin][end].elements
    isempty(prev_el) && return els
    prev_ns = [el[end].numbers for el in els]
    i = 1
    @inbounds while i < lastindex(els[begin])
        for (ie, e) in enumerate(prev_el) 
            j = findfirst(==(e), els[begin][i].elements)
            for (el, prev_n) in zip(els, prev_ns)
                el[i].numbers[j] = el[i].numbers[j] - prev_n[ie]
            end
        end
        i += 1
    end
    els
end