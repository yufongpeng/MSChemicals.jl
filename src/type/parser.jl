abstract type AbstractChemicalParser end
abstract type AbstractAdductParser end

"""
    AdductParser{T} <: AbstractAdductParser

Default adduct parser. `T` determines whether returning a tuple (`true`) or a named tuple (`false`). Namedtuple is useful for user input from a data table and ionization using `ionize`; tuple can be used for fixed adduct ion constructor interface.

# Fields 
* `chemicalparser::T`: chemical parser for each chemical scheme.

# Constructors
    AdductParser(T::Bool = false, parser::AbstractChemicalParser = ChemicalSchemeParser())
"""
struct AdductParser{T} <: AbstractAdductParser
    chemicalparser::AbstractChemicalParser
end
AdductParser() = AdductParser{false}(ChemicalSchemeParser())
AdductParser(x::Bool) = AdductParser{x}(ChemicalSchemeParser())
AdductParser(parser::AbstractChemicalParser) = AdductParser{false}(parser)
AdductParser(x::Bool, parser::AbstractChemicalParser) = AdductParser{x}(parser)

"""
    FormulaChemicalParser <: AbstractChemicalParser

Default chemical parser which parses input string (chemical formula) into `FormulaChemical`. 

# Fields 
* `property::Vector{Pair{Symbol, Any}}`: additional attributes. These can be mutated, added, or deleted with keyword arguments of `parse_chemical`. Values of `nothing` are deleted.

# Constructors
    FormulaChemicalParser(; kwargs...)
    
`kwargs` are collected into field `property`.
"""
struct FormulaChemicalParser <: AbstractChemicalParser
    property::Vector{Pair{Symbol, Any}}
end
FormulaChemicalParser(; kwargs...) = FormulaChemicalParser(collect(kwargs))

"""
    ChemicalParser <: AbstractChemicalParser

Chemical parser which parses input string (chemical name) into `Chemical`. 

# Fields 
* `property::Vector{Pair{Symbol, Any}}`: additional attributes. It must include `:formula` or `:elements`. These can be mutated, added, or deleted with keyword arguments of `parse_chemical`. Values of `nothing` are deleted.

# Constructors
    ChemicalParser(; kwargs...)

`kwargs` are collected into field `property`.
"""
struct ChemicalParser <: AbstractChemicalParser
    property::Vector{Pair{Symbol, Any}}
end
ChemicalParser(; kwargs...) = ChemicalParser(collect(kwargs))

"""
    ChemicalTransitionParser{T} <: AbstractChemicalParser

Default chemical parser which parses input string into `ChemicalTransition`. 

The input string is regarded as series of chemicals separated by " -> ". Pairs and vectors can also be parsed into `ChemicalTransition`. 

If the input can only be parsed into single chemical entity, it is returned directly without wrapping into `ChemicalTransition`.

# Fields 
* `chemicalparser::T`: chemical parser for each chemical in the transition.

# Constructors
    ChemicalTransitionParser(chemicalparser = ChemicalExpressionParser())
"""
struct ChemicalTransitionParser{T<:AbstractChemicalParser} <: AbstractChemicalParser
    chemicalparser::T
end 
ChemicalTransitionParser() = ChemicalTransitionParser(ChemicalExpressionParser())

"""
    ChemicalExpressionParser{T} <: AbstractChemicalParser

Default chemical parser which parses input string into `AbstractChemicalsSchema`. 

# Fields
* `chemicalparser::T`: chemical parser for core chemical. 
* `charge::Int`: default chemical entity charge.
* `gain::Int`: default chemical gain charge.
* `loss::Int`: default chemical loss charge.
* `entity::Bool`: whether parse input into a chemical entity.
* `scheme::Bool`: whether parse input into a chemical scheme.

The input is parsed recursively.
First check `entity` and `scheme`, if `scheme` and `enetity` are true, run step 1 ~ 8; if only `scheme` is true, run only step 1 ~ 6; if only `entity` is true, run only step 7 ~ 8. 
1. Match keys of `scheme_name()` -> `AbstractChemicalScheme`
2. `"+core` -> `"+[core]\$(gain)"` -> step 4
3. `"-core` -> `"-[core]\$(loss)"` -> step 5
4. `"+[core+scheme...]n+"`, `"+[core+scheme...]n-"`, `"+[core+scheme...]"` -> step 8 -> `ChemicalGain` of a chemical entity.
5. `"-[core+scheme...]n+"`, `"-[core+scheme...]n-"`, `"-[core+scheme...]"` -> step 8 -> `ChemicalLoss` of a chemical entity.
6. `"[scheme...]n+"`, `"[scheme...]n-"`, `"[scheme...]"` -> parse each scheme
7. `core` -> `"[core]\$(charge)` -> step 8
8. `"[core+scheme...]n+"`, `"[core+scheme...]n-"`, `"[core+scheme...]"` -> parse `core` with `chemicalparser`; parse the rest part recursively -> a chemical entity

# Constructors
    ChemicalExpressionParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0)
"""
struct ChemicalExpressionParser{T<:AbstractChemicalParser} <: AbstractChemicalParser
    chemicalparser::T 
    charge::Int 
    gain::Int 
    loss::Int
    entity::Bool 
    scheme::Bool
end
ChemicalExpressionParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0) = ChemicalExpressionParser(chemicalparser, charge, gain, loss, true, true)
"""
    ChemicalEntityParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0) 

A `ChemicalExpressionParser` parsing only chemical entity.
"""
ChemicalEntityParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0) = ChemicalExpressionParser(chemicalparser, charge, gain, loss, true, false)
"""
    ChemicalSchemeParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0) 

A `ChemicalExpressionParser` parsing only chemical scheme.
"""
ChemicalSchemeParser(chemicalparser = FormulaChemicalParser(); charge = 0, gain = 0, loss = 0) = ChemicalExpressionParser(chemicalparser, charge, gain, loss, false, true)