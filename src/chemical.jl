"""
    match_chemical(exp, lib; colexp = :Chemical, collib = :Chemical, fnexp = detectedchemical) -> Table

Match chemicals in `exp` (a `Table` or `Vector`) to chemicals in `lib` (a `Table` or `Vector`). 
The resulting table is `exp` with matched index (column `LibID`), matched chemicals (column `Match`) and other information from `lib`.

The exact chemicals being matched are converted from `exp` using `fnexp`.
"""
function match_chemical(exp, lib; colexp = :Chemical, collib = :Chemical, fnexp = detectedchemical)
    del = Int[]
    libid = Int[]
    exp = hasproperty(exp, colexp) ? exp : Table(; Chemical = exp)
    chemical_exp = fnexp.(getproperty(exp, colexp))
    chemical_lib = hasproperty(lib, collib) ? getproperty(lib, collib) : lib
    for i in eachindex(exp)
        j = findfirst(x -> ischemicalequal(x, chemical_exp[i]), chemical_lib)
        isnothing(j) ? push!(del, i) : push!(libid, j)
    end
    id = setdiff(eachindex(exp), del)
    ps = filter(!=(collib), propertynames(lib))
    Table(exp[id]; LibID = libid, Match = chemical_lib[libid], [p => [getproperty(lib, p)[i] for i in libid] for p in ps]...)
end

"""
    ischemicalequal(x::AbstractChemicalsSchema, y::AbstractChemicalsSchema) -> Bool

Determine whether two chemicals are chemically equivalent. 
By default, it transforms both chemicals by [`ischemicalequaltransform`](@ref) and compares them by [`istransformedchemicalequal`](@ref).
"""
ischemicalequal(x::AbstractChemical, y::AbstractChemical) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))
ischemicalequal(x::AbstractScheme, y::AbstractScheme) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))
ischemicalequal(x::Isobars, y::Isobars) = istransformedchemicalequal(x, y)
ischemicalequal(x::Isotopomers, y::Isotopomers) = istransformedchemicalequal(x, y)
ischemicalequal(x::Groupedisotopomers, y::Groupedisotopomers) = istransformedchemicalequal(x, y)
ischemicalequal(x::ChemicalTransition, y::ChemicalTransition) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))
ischemicalequal(x::ChemicalSchema, y::ChemicalSchema) = istransformedchemicalequal(x, y)
ischemicalequal(x::StructuralElementalScheme, y::StructuralElementalScheme) = istransformedchemicalequal(x, y)
ischemicalequal(x::ElementalScheme{true}, y::ElementalScheme{true}) = istransformedchemicalequal(x, y)
ischemicalequal(x::ElementalScheme{false}, y::ElementalScheme{false}) = istransformedchemicalequal(x, y)
ischemicalequal(x::IsotopomerizedSchema, y::IsotopomerizedSchema) = istransformedchemicalequal(x, y)
ischemicalequal(x::Groupedisotopomerizedschema, y::Groupedisotopomerizedschema) = istransformedchemicalequal(x, y)
ischemicalequal(x::AbstractChemicalsSchema, y::AbstractChemicalsSchema) = istransformedchemicalequal(ischemicalequaltransform(x), ischemicalequaltransform(y))

"""
    ischemicalequaltransform(x::AbstractChemicalsSchema) -> AbstractChemicalsSchema

Return an object for comparison with other chemicals by [`istransformedchemicalequal`](@ref). 
"""
ischemicalequaltransform(x::AbstractChemical) = x 
ischemicalequaltransform(x::AbstractScheme) = x 
ischemicalequaltransform(x::T) where {T<:AbstractChemicalWrapper} = ischemicalequaltransform(x.chemical)
ischemicalequaltransform(x::Isobars) = length(x) == 1 ? ischemicalequaltransform(chemicalentity(x)) : x
ischemicalequaltransform(x::Isotopomers) = isempty(unique_elements(x.isotopes)) ? x.parent : x 
ischemicalequaltransform(x::Groupedisotopomers) = length(x.isotopes) > 1 ? x : isempty(unique_elements(x.isotopes[begin])) ? x.parent : Isotopomers(x.parent, x.isotopes[begin]) 
ischemicalequaltransform(x::IsotopomerizedSchema) = isempty(unique_elements(x.isotopes)) ? x.parent : x 
ischemicalequaltransform(x::Groupedisotopomerizedschema) = length(x.isotopes) > 1 ? x : isempty(unique_elements(x.isotopes[begin])) ? x.parent : IsotopomerizedSchema(x.parent, x.isotopes[begin]) 
ischemicalequaltransform(x::ElementalScheme{T}) where T = ElementalScheme(T, ischemicalequaltransform(x.chemical))
ischemicalequaltransform(x::StructuralElementalScheme) = StructuralElementalScheme(ischemicalequaltransform(x.structuralscheme), ischemicalequaltransform(x.elementalscheme))
ischemicalequaltransform(x::ChemicalTransition) = ChemicalTransition([ischemicalequaltransform(c) for c in chemicaltransition(x)])

"""
    istransformedchemicalequal(x::AbstractChemicalsSchema, y::AbstractChemicalsSchema) -> Bool

Determine whether two chemicals are chemically equivalent after applying [`ischemicalequaltransform`](@ref). 
For [`Chemical`](@ref) and [`FormulaChemical`](@ref), It tests the name and the elements composition.
"""
istransformedchemicalequal(x::AbstractChemicalsSchema, y::AbstractChemicalsSchema) = false
istransformedchemicalequal(x::AbstractChemical, y::AbstractChemical) = isequal(x, y)
istransformedchemicalequal(x::AbstractScheme, y::AbstractScheme) = isequal(x, y)
istransformedchemicalequal(x::AbstractAdductIon, y::AbstractAdductIon) = ischemicalequal(ionadduct(x), ionadduct(y)) && ischemicalequal(ioncore(x), ioncore(y))
istransformedchemicalequal(x::Chemical, y::Chemical) = 
    isequal(chemicalname(x), chemicalname(y)) && isequal(sort_unique_elements(chemicalelements(x)), sort_unique_elements(chemicalelements(y)))
istransformedchemicalequal(x::FormulaChemical, y::FormulaChemical) = 
    isequal(chemicalname(x), chemicalname(y)) 
istransformedchemicalequal(x::Isobars, y::Isobars) = all(ischemicalequal(a, b) for (a, b) in zip(x.chemicals, y.chemicals)) && all(isapprox(a, b) for (a, b) in zip(x.abundance, y.abundance))
istransformedchemicalequal(x::Isotopomers, y::Isotopomers) = ischemicalequal(x.parent, y.parent) && x.isotopes == y.isotopes
istransformedchemicalequal(x::Groupedisotopomers, y::Groupedisotopomers) = ischemicalequal(x.parent, y.parent) && x.state == y.state && x.isotope == y.isotope && all(splat(==), zip(x.isotopes, y.isotopes)) && all(splat(isapprox), zip(x.abundance, y.abundance))
istransformedchemicalequal(x::ChemicalTransition, y::ChemicalTransition) = all(ischemicalequal.(x.transition, y.transition))
istransformedchemicalequal(x::IsotopomerizedSchema, y::IsotopomerizedSchema) = istransformedchemicalequal(x.parent, y.parent) && x.isotopes == y.isotopes
istransformedchemicalequal(x::Groupedisotopomerizedschema, y::Groupedisotopomerizedschema) = ischemicalequal(x.parent, y.parent) && x.state == y.state && x.isotope == y.isotope && all(splat(==), zip(x.isotopes, y.isotopes)) && all(splat(isapprox), zip(x.abundance, y.abundance))
function istransformedchemicalequal(x::ChemicalSchema, y::ChemicalSchema) 
    uk = [false for _ in eachindex(y.schema)]
    for (kx, vx) in zip(x.schema, x.number)
        pass = false
        for (i, ky) in enumerate(y.schema)
            uk[i] && continue 
            if ischemicalequal(kx, ky) && vx == y.number[i]
                pass = true
                uk[i] = true
                break
            end
        end
        pass || return false
    end
    all(uk)
end
istransformedchemicalequal(x::ElementalScheme, y::ElementalScheme) = ischemicalequal(x.chemical, y.chemical)
istransformedchemicalequal(x::StructuralElementalScheme, y::StructuralElementalScheme) = ischemicalequal(structuralscheme(x), structuralscheme(y)) && ischemicalequal(elementalscheme(x), elementalscheme(y))

"""
    ionize([constructor,] chemical; kwargs...) -> constructor
    ionize([constructor,] chemical, args...; kwargs...) -> constructor

Ionize `chemical` and wrap with `constructor`. The default constructor is [`AdductIon`](@ref).
"""
ionize(chemical::AbstractChemical; kwargs...) = ionize(AdductIon, chemical; kwargs...)
ionize(chemical::AbstractChemical, adduct, ncore = 1; kwargs...) = ionize(AdductIon, chemical, adduct, ncore; kwargs...)

"""
    ionize(::Type{AdductIon}, chemical::AbstractChemical; adduct, ncore = 1, kwargs...) -> AdductIon
    ionize(::Type{AdductIon}, chemical::AbstractChemical, adduct, ncore = 1; kwargs...) -> AdductIon

Ionize `chemical` and wrap with [`AdductIon`](@ref).

# Arguments 
* `adduct::AbstractScheme`: adduct scheme.
* `ncore::Int`: number of core chemicals.
"""
ionize(::Type{AdductIon}, chemical::AbstractChemical; adduct, ncore = 1, kwargs...) = AdductIon(chemical, adduct, ncore)
ionize(::Type{AdductIon}, chemical::AbstractAdductIon; adduct, ncore = 1, kwargs...) = AdductIon(ioncore(chemical), adduct, ncore)
ionize(::Type{AdductIon}, chemical::Isobars; adduct, ncore = 1, kwargs...) = Isobars([ionize(AdductIon, x, adduct, ncore; kwargs...) for x in chemicalspecies(chemical)], chemical.abundance)
ionize(::Type{AdductIon}, chemical::Isotopomers; adduct, ncore = 1, kwargs...) = Isotopomers(ionize(AdductIon, chemicalparent(chemical), adduct, ncore; kwargs...), chemical.isotopes)
ionize(::Type{AdductIon}, chemical::Groupedisotopomers; adduct, ncore = 1, kwargs...) = Groupedisotopomers(ionize(AdductIon, chemicalparent(chemical), adduct, ncore; kwargs...), chemical.state, chemical.isotope, chemical.isotopes, chemical.abundance)

ionize(::Type{AdductIon}, chemical::AbstractChemical, adduct, ncore = 1; kwargs...) = AdductIon(chemical, adduct, ncore)
ionize(::Type{AdductIon}, chemical::AbstractAdductIon, adduct, ncore = 1; kwargs...) = AdductIon(ioncore(chemical), adduct, ncore)
ionize(::Type{AdductIon}, chemical::Isobars, adduct, ncore = 1; kwargs...) = Isobars([ionize(AdductIon, x, adduct, ncore; kwargs...) for x in chemicalspecies(chemical)], chemical.abundance)
ionize(::Type{AdductIon}, chemical::Isotopomers, adduct, ncore = 1; kwargs...) = Isotopomers(ionize(AdductIon, chemicalparent(chemical), adduct, ncore; kwargs...), chemical.isotopes)
ionize(::Type{AdductIon}, chemical::Groupedisotopomers, adduct, ncore = 1; kwargs...) = Groupedisotopomers(ionize(AdductIon, chemicalparent(chemical), adduct, ncore; kwargs...), chemical.state, chemical.isotope, chemical.isotopes, chemical.abundance)

"""
    isotopomerize(chemical::AbstractChemicalsSchema, isotopes) -> AbstractChemicalsSchema

Add delocalized isotopic replacements `isotopes` to `chemical`.
"""
isotopomerize(chemical::AbstractChemical, isotopes) = Isotopomers(chemical, isotopes)
isotopomerize(chemical::Isotopomers, isotopes) = Isotopomers(chemical.parent, gain_elements(chemical.isotopes, isotopes))
isotopomerize(sch::StructuralElementalScheme, isotopes) = StructuralElementalScheme(structuralscheme(sch), isotopomerize(elementalscheme(sch), isotopes))
# isotopomerize(sch::AbstractCompleteScheme{T,<:AbstractChemical}, isotopes) where T = isotopomerize(elementalscheme(sch), isotopes)
# isotopomerize(sch::StructuralElementalScheme{T,<:AbstractChemical}, isotopes) where T = isotopomerize(elementalscheme(sch), isotopes)
isotopomerize(sch::ElementalScheme{true}, isotopes) = ElementalScheme(true, isotopomerize(sch.chemical, isotopes))
isotopomerize(sch::ElementalScheme{false}, isotopes) = ElementalScheme(false, isotopomerize(sch.chemical, reverse_elements(isotopes, true)))
# isotopomerize(sch::ElementalScheme{false}, isotopes::Tuple) = ElementalScheme(false, isotopomerize(sch.chemical, loss_elements(isotopes...)))
isotopomerize(sch::ChemicalSchema, isotopes) = IsotopomerizedSchema(sch, isotopes)
# isotopomerize(sch::ChemicalSchema, isotopes::Tuple) = IsotopomerizedSchema(sch, loss_elements(last(isotopes), first(isotopes)))
isotopomerize(sch::IsotopomerizedSchema, isotopes) = IsotopomerizedSchema(sch.parent, gain_elements(sch.isotopes, isotopes))
# isotopomerize(sch::IsotopomerizedSchema, isotopes::Tuple) = IsotopomerizedSchema(sch.parent, loss_elements!(gain_elements(sch.isotopes, last(isotopes), first(isotopes))))
isotopomerize(sch::T, isotopes) where {T<:AbstractScheme} = throw(ArgumentError("Cannot add isotopes information to $T."))

"""
    groupedisotopomerize(chemical::AbstractChemicalsSchema, isotopes) -> AbstractChemicalsSchema

Add delocalized isotopic replacements `isotopes` to `chemical`.
"""
groupedisotopomerize(chemical::AbstractChemical, state, isotope, isotopes, abundance) = Groupedisotopomers(chemical, state, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(chemical::Isotopomers, state, isotope, isotopes, abundance) = Groupedisotopomers(chemical.parent, state, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(chemical::Groupedisotopomers, state, isotope, isotopes, abundance) = Groupedisotopomers(chemical.parent, state, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(sch::StructuralElementalScheme, state, isotope, isotopes, abundance) = StructuralElementalScheme(structuralscheme(sch), groupedisotopomerize(elementalscheme(sch), state, isotope, isotopes, abundance))
# groupedisotopomerize(sch::CompleteSchemeChemical, state, isotope, isotopes, abundance) = groupedisotopomerize(elementalscheme(sch), state, isotope, isotopes, abundance)
# groupedisotopomerize(sch::StructuralElementalScheme{T,<:AbstractChemical}, state, isotope, isotopes, abundance) where T = groupedisotopomerize(elementalscheme(sch), state, isotope, isotopes, abundance)
groupedisotopomerize(sch::ElementalScheme{true}, state, isotope, isotopes, abundance) = ElementalScheme(true, groupedisotopomerize(sch.chemical, state, isotope, isotopes, abundance))
groupedisotopomerize(sch::ElementalScheme{false}, state, isotope, isotopes, abundance) = ElementalScheme(false, groupedisotopomerize(sch.chemical, state, isotope, reverse_elements.(isotopes, true), abundance))
groupedisotopomerize(sch::ChemicalSchema, state, isotope, isotopes, abundance) = Groupedisotopomerizedschema(sch, state, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(sch::IsotopomerizedSchema, state, isotope, isotopes, abundance) = Groupedisotopomerizedschema(sch.parent, state, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(sch::Groupedisotopomerizedschema, state, isotope, isotopes, abundance) = Groupedisotopomerizedschema(sch.parent, state, isotope, groupedisotopomersisotopes(isotopes), abundance)
groupedisotopomerize(sch::T, state, isotope, isotopes, abundance) where {T<:AbstractScheme} = throw(ArgumentError("Cannot add isotopes information to $T."))

groupedisotopomersisotopes(isotopes::Vector{ElementsVector}) = isotopes
groupedisotopomersisotopes(isotopes::Vector{Vector{Pair{String, Int64}}}) = [ElementsVector(first.(x), last.(x)) for x in isotopes]
groupedisotopomersisotopes(isotopes::Vector{Dict{String, Int64}}) = [ElementsVector(collect(keys(x)), collect(values((x)))) for x in isotopes]
groupedisotopomersisotopes(isotopes::Vector{Dictionary{String, Int64}}) = [ElementsVector(collect(keys(x)), collect(values((x)))) for x in isotopes]

"""
    isgainscheme(sch::AbstractScheme) -> Bool

Whether `sch` contains only chemical gains.
"""
isgainscheme(sch) = false
isgainscheme(sch::ElementalScheme{true}) = true
isgainscheme(sch::AbstractCompleteScheme) = isgainscheme(elementalscheme(sch))
isgainscheme(sch::ChemicalSchema) = all(isgainscheme, sch.schema)
isgainscheme(sch::IsotopomerizedSchema) = isgainscheme(sch.parent)
isgainscheme(sch::Groupedisotopomerizedschema) = isgainscheme(sch.parent)

"""
    islossscheme(sch::AbstractScheme) -> Bool

Whether `sch` contains only chemical losses.
"""
islossscheme(sch) = false
islossscheme(sch::ElementalScheme{false}) = true
islossscheme(sch::AbstractCompleteScheme) = islossscheme(elementalscheme(sch))
islossscheme(sch::ChemicalSchema) = all(islossscheme, sch.schema)
islossscheme(sch::IsotopomerizedSchema) = islossscheme(sch.parent) 
islossscheme(sch::Groupedisotopomerizedschema) = islossscheme(sch.parent) 

"""
    completescheme(precursor::AbstractChemical, product::AbstractChemical) -> StructuralElementalScheme
    completescheme(precursor::AbstractChemical, sch::AbstractScheme) -> StructuralElementalScheme    
    completescheme(precursor::AbstractChemical, product::GenericChemical) -> StructuralElementalScheme
    completescheme(precursor::AbstractChemical, product::AdductIon{<:GenericChemical}) -> StructuralElementalScheme
    completescheme(precursor::AbstractChemical, sch::CompleteSchema) -> CompleteSchema
    completescheme(precursor::AbstractChemical, sch::ChemicalSchema) -> ChemicalSchema 
    completescheme(precursor::AbstractChemical, sch::IsotopomerizedSchema) -> IsotopomerizedSchema
    completescheme(precursor::AbstractChemical, sch::Groupedisotopomerizedschema) -> Groupedisotopomerizedschema
    completescheme(precursor::Nothing, product::AbstractChemical) -> StructuralElementalScheme

Transform `sch` or `product` into `CompleteSchema` according to `precursor`. 

For generic precursor and product, this function calls `structure_search` which searches the property `:structure` of `ioncore(precursor)` for `ionadduct(precursor)` and `product`. 
When `product` is a generic chemical, `structure_search` searches the property `:chemicalscheme` for `ionadduct(precursor)` first, and then use the returned scheme instead of `product` for the following structure search. 

For other chemicals, it calls `elementalscheme(precursor, product)` to generate elemental scheme incorporating information from `precursor`. 

# Examples
## Generic chemical
```julia-repl
julia> chemical = Chemical("18:0 PC-d9", "C44H79NO8PD9") # Deuterium-labeled PC on methyl group
18:0 PC-d9

julia> loss_me = ElementalScheme(false, Chemical("Me", "CH3"; charge = 1))
Loss_Me

julia> loss_cd3 = ElementalScheme(false, Chemical("Me[D3]", "CD3"; charge = 1))
Loss_Me[D3]

julia> push!(chemical.property, :structure => [nothing => [loss_me => loss_cd3]])
1-element Vector{Pair{Symbol, Any}}:
 :structure => Pair{Nothing, Vector{Pair{ElementalScheme{false, Chemical}, ElementalScheme{false, Chemical}}}}[nothing => [Loss_Me => Loss_Me[D3]]]

julia> ionize(chemical, loss_me)
[(18:0 PC-d9)-Me[D3]]-
```
A more complex example, 
```julia-repl
julia> ps1 = Chemical("PS[D3,13C3] 18:0/20:4", "C41[13C]3H75D3NO10P") # Deuterium/carbon-13-labeled PS on serine part
PS[D3,13C3] 18:0/20:4

julia> ps2 = Chemical("PS 18:0[D5]/20:4(5Z,8Z,11Z,14Z)", "C44H73D5NO10P") # Deuterium/carbon-13-labeled PS on fa1 part
PS 18:0[D5]/20:4(5Z,8Z,11Z,14Z)

julia> begin
       cserine = Chemical("Serine", "C3H5NO2"; abbreviation = "Ser") # Normal serine
       cserinei = Chemical("Serine[D3,13C3]", "[13C]3H2D3NO2"; abbreviation = "Ser[D3,13C3]") # Labeled serine
       lossserine = ChemicalLoss(cserine)
       lossserinei = ChemicalLoss(cserinei)
       losshserine = ChemicalLoss(AdductIon(cserine, "[M+H]+"))
       losshserinei = ChemicalLoss(AdductIon(cserinei, "[M+H]+"))
       fa1 = Chemical("FA 18:0", "C18H36O2")
       fa1i = Chemical("FA 18:0[D5]", "C13D5H36O2") # Labeled fa
       end;

julia> push!(ps1.property, :structure => [
           nothing => [
               lossserine => lossserinei,
               losshserine => losshserinei
           ],
           ChemicalLoss(Proton()) => [
               lossserine => lossserinei  
           ]
       ]);

julia> push!(ps1.property, :schema => [
           ChemicalLoss(Proton()) => [
               lossserine => losshserine
           ] 
       ]);

julia> detectedchemical(ChemicalSeries(ionize(ps1, ChemicalLoss(Proton())) => lossserine))
[(PS[D3,13C3] 18:0/20:4)-Ser[D3,13C3]-H]-

julia> push!(fa1.property, :chemicalscheme => [
           ChemicalLoss(Proton()) => :sn1fa
       ]);

julia> push!(fa1i.property, :chemicalscheme => [
           ChemicalLoss(Proton()) => :sn1fa
       ]);

julia> push!(ps2.property, :structure => [
           ChemicalLoss(Proton()) => [
               :sn1fa => AdductIon(fa1i, ChemicalLoss(Proton()))     
           ],
           losshserine => [
               :sn1fa => AdductIon(fa1i, ChemicalLoss(Proton()))     
           ]
       ]);
    
julia> push!(ps2.property, :schema => [
           ChemicalLoss(Proton()) => [
               lossserine => losshserine
           ] 
       ]);

julia> detectedchemical(ChemicalSeries(ionize(ps2, ChemicalLoss(Proton())) => fa1))
[FA 18:0[D5]-H]-
```
## Customized type
```julia-repl
julia> abstract type AbstractPC <: AbstractChemical end

julia> struct PC <: AbstractPC end # Normal PC

julia> struct DLPC <: AbstractPC 
           location::Symbol
       end # Deuterium-labeled PC on methyl group (location = :Me) or other part

julia> struct Me <: AbstractChemical end; charge(::Me; kwargs...) = -1.0 # Methenium
charge (generic function with 15 methods)

julia> struct DLMe <: AbstractChemical end; charge(::DLMe; kwargs...) = -1.0 # Deuterium-labeled Methinium
charge (generic function with 16 methods)

julia> elementalscheme(pc::DLPC, ::CLType(Me)) = pc.location == :Me ? ChemicalLoss(DLMe()) : ChemicalLoss(Me())
elementalscheme (generic function with 31 methods)

julia> ion1 = ionize(DLPC(:Me); adduct = ChemicalLoss(Me()))
[DLPC-DLMe]-

julia> adductionscheme(pc::AIType(AbstractPC, SESType(CGType(Acetate))), ::CLType(MethylAcetate)) = completescheme(ioncore(pc), ChemicalLoss(Me()))
adductionscheme (generic function with 5 methods)

julia> ion2 = detectedchemical(ChemicalSeries(ionize(DLPC(:Me), "[M+OAc]-") => ChemicalLoss(MethylAcetate())))
[DLPC-DLMe]-

julia> ion1 == ion2
true
```
"""
completescheme(precursor::AbstractChemical, product::AbstractChemical) = StructuralElementalScheme(RandomProductScheme(), product)
completescheme(precursor::Nothing, product::AbstractChemical) = StructuralElementalScheme(RandomProductScheme(), product)
# completescheme(precursor::T, product::S) where {T<:AbstractChemical, S<:AbstractScheme} = throw(ArgumentError("Specific `completescheme(precursor::$T, scheme::$S)` method has to be implemented."))
completescheme(precursor::AbstractChemical, product::CompleteSchema) = product
completescheme(precursor::AbstractChemical, product::GenericChemical) = _completescheme(precursor, product)
completescheme(precursor::AbstractChemical, product::AdductIon{<:GenericChemical}) = _completescheme(precursor, product)
completescheme(precursor::AbstractChemical, product::ChemicalSchema) = _completescheme(precursor, product)
completescheme(precursor::AbstractChemical, product::IsotopomerizedSchema) = _completescheme(precursor, product)
completescheme(precursor::AbstractChemical, product::Groupedisotopomerizedschema) = _completescheme(precursor, product)
# completescheme(precursor::AbstractChemical, product::AbstractElementalScheme) = StructuralElementalScheme(product, copy(product))
completescheme(precursor::AbstractChemical, product::AbstractScheme) = StructuralElementalScheme(product, elementalscheme(precursor, product))

# _completescheme(precursor::AbstractChemical, product::ElementalScheme{T}) where T = StructuralElementalScheme(product, ElementalScheme(T, product))

_completescheme(precursor::AbstractChemical, product::GenericChemical) = StructuralElementalScheme(RandomProductScheme(), product)
_completescheme(precursor::AbstractChemical, product::AdductIon{<:GenericChemical}) = StructuralElementalScheme(RandomProductScheme(), product)
_completescheme(precursor::AbstractChemical, product::ChemicalSchema) = ChemicalSchema(completescheme.(Ref(precursor), product.schema), product.number)
_completescheme(precursor::AbstractChemical, product::IsotopomerizedSchema) = IsotopomerizedSchema(completescheme(precursor, product.parent), product.isotopes)
_completescheme(precursor::AbstractChemical, product::Groupedisotopomerizedschema) = Groupedisotopomerizedschema(completescheme(precursor, product.parent), product.state, product.isotope, product.isotopes, product.abundance)
# structure search for generic types
_completescheme(precursor::GenericChemical, product::GenericChemical) = structure_search(precursor, nothing, product)
_completescheme(precursor::GenericChemical, product::AdductIon{<:GenericChemical}) = structure_search(precursor, nothing, product)
# _completescheme(precursor::GenericChemical, product::AbstractScheme) = structure_search(precursor, nothing, product)
# _completescheme(precursor::GenericChemical, product::CompleteSchema) = product
_completescheme(precursor::GenericChemical, product::ChemicalSchema) = structure_search(precursor, nothing, product)
_completescheme(precursor::GenericChemical, product::IsotopomerizedSchema) = structure_search(precursor, nothing, product)
_completescheme(precursor::GenericChemical, product::Groupedisotopomerizedschema) = structure_search(precursor, nothing, product)
# _completescheme(precursor::GenericChemical, product::AbstractElementalScheme) = structure_search(precursor, nothing, product)
_completescheme(precursor::AdductIon{<:GenericChemical}, product::GenericChemical) = structure_search(ioncore(precursor), ionadduct(precursor), product)
_completescheme(precursor::AdductIon{<:GenericChemical}, product::AdductIon{<:GenericChemical}) = structure_search(ioncore(precursor), ionadduct(precursor), product)
# _completescheme(precursor::AdductIon{<:GenericChemical}, product::AbstractScheme) = structure_search(ioncore(precursor), ionadduct(precursor), product)
# _completescheme(precursor::AdductIon{<:GenericChemical}, product::CompleteSchema) = product
_completescheme(precursor::AdductIon{<:GenericChemical}, product::IsotopomerizedSchema) = structure_search(ioncore(precursor), ionadduct(precursor), product)
_completescheme(precursor::AdductIon{<:GenericChemical}, product::Groupedisotopomerizedschema) = structure_search(ioncore(precursor), ionadduct(precursor), product)
_completescheme(precursor::AdductIon{<:GenericChemical}, product::ChemicalSchema) = structure_search(ioncore(precursor), ionadduct(precursor), product)
# _completescheme(precursor::AdductIon{<:GenericChemical}, product::AbstractElementalScheme) = structure_search(ioncore(precursor), ionadduct(precursor), product)

"""
    completeschemechemical(precursor::AbstractChemicalsSchema, product::AbstractChemicalsSchema) -> AbstractChemicalsSchema
    completeschemechemical(sch::CompleteSchemeChemical) -> AbstractChemical
    completeschemechemical(sch::AbstractScheme) -> AbstractScheme

Transform `product` into `CompleteSchema` or `AbstractChemical` according to `precursor`. `[`completescheme`](@ref)` is first called to generate `sch`, and the returned object is determined by the type of `sch`. If `sch` is a `CompleteSchemeChemial`, then `elementalscheme(sch)` is called; otherwise, `sch` is returned as is.
"""
completeschemechemical(precursor, product) = completeschemechemical(completescheme(precursor, product))
completeschemechemical(sch::CompleteSchemeChemical) = elementalscheme(sch)
completeschemechemical(sch::AbstractScheme) = sch

"""
    elementalscheme(precursor::AbstractChemical, product::AbstractChemical) -> AbstractChemical
    elementalscheme(precursor::AbstractChemical, sch::AbstractScheme) -> StructuralElementalScheme    
    elementalscheme(precursor::AbstractChemical, sch::CompleteSchema) -> CompleteSchema
    elementalscheme(precursor::AbstractChemical, sch::ChemicalSchema) -> ChemicalSchema 
    elementalscheme(precursor::AbstractChemical, sch::IsotopomerizedSchema) -> IsotopomerizedSchema
    elementalscheme(precursor::AbstractChemical, sch::Groupedisotopomerizedschema) -> Groupedisotopomerizedschema
    elementalscheme(precursor::AbstractChemical, sch::AbstractElementalScheme) -> AbstractElementalScheme
    elementalscheme(precursor::Nothing, product::AbstractChemical) -> AbstractChemical

Transform `sch` or `product` into elemental scheme according to `precursor`. 

For generic precursor and product, this function calls `structure_search_elemental` which searches the property `:structure` of `ioncore(precursor)` for `ionadduct(precursor)` and `product`. 
When `product` is a generic chemical, `structure_search_elemental` searches the property `:chemicalscheme` for `ionadduct(precursor)` first, and then use the returned scheme instead of `product` for the following structure search. 

Any of the following method should be defined for new structural scheme and new chemical type:
* `elementalscheme(::new_chemical_type, ::new_structural_type)`: the new chemical type represents an ion in MS.
* `elementalscheme(::AdductIon{new_chemical_type, StructuralElementalScheme{structural_type}}, ::new_structural_type)`: the new chemical type formed an adduct ion with single scheme.
* `elementalscheme(::AdductIon{new_chemical_type, ChemicalSchema}, ::new_structural_type)`: the new chemical type formed an adduct ion with multiple schema.

See [`completescheme`](@ref) for examples of defining new methods.
```
"""
elementalscheme(precursor::AbstractChemical, product::AbstractChemical) = product
elementalscheme(precursor::Nothing, product::AbstractChemical) = product
elementalscheme(precursor::T, product::S) where {T<:AbstractChemical, S<:AbstractScheme} = throw(ArgumentError("Specific `elementalscheme(precursor::$T, scheme::$S)` method has to be implemented."))
elementalscheme(precursor::AbstractChemical, product::CompleteSchema) = elementalscheme(product)
elementalscheme(precursor::AbstractChemical, product::ChemicalSchema) = ChemicalSchema(elementalscheme.(Ref(precursor), product.schema), product.number)
elementalscheme(precursor::AbstractChemical, product::IsotopomerizedSchema) = IsotopomerizedSchema(elementalscheme(precursor, product.parent), product.isotopes)
elementalscheme(precursor::AbstractChemical, product::Groupedisotopomerizedschema) = Groupedisotopomerizedschema(elementalscheme(precursor, product.parent), product.state, product.isotope, product.isotopes, product.abundance)
elementalscheme(precursor::AbstractChemical, product::AbstractElementalScheme) = copy(product)
# elementalscheme(precursor::AbstractChemical, product::ElementalScheme{T}) where T = StructuralElementalScheme(product, ElementalScheme(T, product))
# structure search for generic types
elementalscheme(precursor::GenericChemical, product::GenericChemical) = structure_search_elemental(precursor, nothing, product)
elementalscheme(precursor::GenericChemical, product::AdductIon{<:GenericChemical}) = structure_search_elemental(precursor, nothing, product)
elementalscheme(precursor::GenericChemical, product::AbstractScheme) = structure_search_elemental(precursor, nothing, product)
elementalscheme(precursor::GenericChemical, product::CompleteSchema) = elementalscheme(product)
elementalscheme(precursor::GenericChemical, product::ChemicalSchema) = structure_search_elemental(precursor, nothing, product)
elementalscheme(precursor::GenericChemical, product::IsotopomerizedSchema) = structure_search_elemental(precursor, nothing, product)
elementalscheme(precursor::GenericChemical, product::Groupedisotopomerizedschema) = structure_search_elemental(precursor, nothing, product)
elementalscheme(precursor::GenericChemical, product::AbstractElementalScheme) = structure_search_elemental(precursor, nothing, product)
elementalscheme(precursor::AdductIon{<:GenericChemical}, product::GenericChemical) = structure_search_elemental(ioncore(precursor), ionadduct(precursor), product)
elementalscheme(precursor::AdductIon{<:GenericChemical}, product::AdductIon{<:GenericChemical}) = structure_search_elemental(ioncore(precursor), ionadduct(precursor), product)
elementalscheme(precursor::AdductIon{<:GenericChemical}, product::AbstractScheme) = structure_search_elemental(ioncore(precursor), ionadduct(precursor), product)
elementalscheme(precursor::AdductIon{<:GenericChemical}, product::CompleteSchema) = elementalscheme(product)
elementalscheme(precursor::AdductIon{<:GenericChemical}, product::ChemicalSchema) = structure_search_elemental(ioncore(precursor), ionadduct(precursor), product)
elementalscheme(precursor::AdductIon{<:GenericChemical}, product::IsotopomerizedSchema) = structure_search_elemental(ioncore(precursor), ionadduct(precursor), product)
elementalscheme(precursor::AdductIon{<:GenericChemical}, product::Groupedisotopomerizedschema) = structure_search_elemental(ioncore(precursor), ionadduct(precursor), product)
elementalscheme(precursor::AdductIon{<:GenericChemical}, product::AbstractElementalScheme) = structure_search_elemental(ioncore(precursor), ionadduct(precursor), product)

"""
    adductionscheme(precursor::AdductIon, product::AbstractScheme) -> CompleteSchema

Return a new scheme blending `ionadduct(precursor)` and `product`. 

For generic chemical types, this function calls `schema_search` which searches the property `:schema` of `ioncore(precursor)` for `ionadduct(precursor)` and `product`, 
and then calls [`structure_search`](@ref), seaching for `nothing` and the returned schema. 

For other chemicals, it returns a complete scheme directly without incorporating any information from `precursor`.

Defining new method is optional for new structural scheme and new chemical type unless schema have to be blended.
* `adductionscheme(::AdductIon{new_chemical_type, StructuralElementalScheme{structural_type}}, ::new_structural_type)`.
* `adductionscheme(::AdductIon{new_chemical_type, ChemicalSchema}, ::new_structural_type)`.
* `adductionscheme(::AdductIon{new_chemical_type, StructuralElementalScheme{structural_type}}, ::ChemicalSchema)`.
* `adductionscheme(::AdductIon{new_chemical_type, ChemicalSchema}, ::ChemicalSchema)`.

See [`completescheme`](@ref) for examples of defining new methods.
```
"""
adductionscheme(precursor::AdductIon, product::CompleteSchema) = adductionscheme(precursor, structuralscheme(product))
adductionscheme(precursor::AdductIon, product::AbstractScheme) = ChemicalSchema(ionadduct(precursor), completescheme(precursor, product))
# adductionscheme(precursor::AdductIon, product::CompleteSchema) = StructuralElementalScheme(ChemicalSchema(structuralscheme(ionadduct(precursor)), structuralscheme(product)), ChemicalSchema(elementalscheme(ionadduct(precursor)), elementalscheme(product)))
# schema search for generic adduction
function adductionscheme(precursor::AdductIon{<:GenericChemical}, product::CompleteSchema) 
    sch = schema_search(ioncore(precursor), ionadduct(precursor), product)
    isnothing(sch) ? ChemicalSchema(ionadduct(precursor), product) : structure_search(ioncore(precursor), nothing, sch)
end
function adductionscheme(precursor::AdductIon{<:GenericChemical}, product::AbstractScheme) 
    sch = schema_search(ioncore(precursor), ionadduct(precursor), product)
    isnothing(sch) ? ChemicalSchema(ionadduct(precursor), completescheme(precursor, product)) : structure_search(ioncore(precursor), nothing, sch)
end

# For customized adduction type of specific chemical type, implement 
# detectedchemical(::chemicaltype, ::CompleteSchema) -> adductiontype
# detectedchemical(::adductiontype, ::CompleteSchema) -> adductiontype
# adductionscheme(::adductiontype, ::CompleteSchema) -> CompleteSchema
# completescheme(::chemicaltype, ::AbstractScheme) -> CompleteSchema
# completescheme(::adductiontype, ::AbstractScheme) -> CompleteSchema

# Internal interfaces for detectedchemical
"""
    detectedchemical(precursor::AbstractChemical, product::AbstractChemical) -> AbstractChemical 
    detectedchemical(precursor::AbstractChemical, sch::AbstractScheme) -> AbstractChemical 
    detectedchemical(precursor::AbstractChemical, sch::CompleteSchemeChemical) -> AbstractChemical 
    detectedchemical(precursor::AbstractChemical, sch::CompleteSchema) -> AbstractChemical 
    detectedchemical(precursor::AbstractChemical, sch::StructuralChemicalScheme) -> AbstractChemical 
    detectedchemical(precursor::Nothing, product::AbstractChemicalsSchema) -> AbstractChemical 

The chemical directly detected in MS. 

By default, [`completescheme`](@ref) is applied to `sch` first, and [`AdductIon`](@ref) is generated; chemical entity of `product` or chemical containing `sch` is directly returned.

When `precursor` is [`Isotopomers`](@ref) or [`Groupedisotopomers`](@ref), the ouput chemical entity is wrapped accordingly. 

For [`AdductIon`](@ref), [`adductionscheme`](@ref) is called for blending precursor and product scheme. 

Defining new method is optional unless other [`AbstractAdductIon`](@ref) type is used.
* `detectedchemical(::new_adduction_type, ::CompleteSchema)`.
* `detectedchemical(::new_adduction_type, ::AbstractScheme)`.
For [`StructuralChemicalScheme`](@ref) method, defining new method `elementalscheme(::new_adduction_type, ::structural_type)` for each `structural_type<:StructuralChemicalScheme`.
"""
detectedchemical(precursor::AbstractChemical, product::AbstractChemical) = detectedchemical(precursor, completescheme(precursor, product))
# detectedchemical(precursor::AbstractChemical, product::AbstractScheme) = detectedchemical(precursor, completescheme(precursor, product))
detectedchemical(precursor::AbstractChemical, product::CompleteSchemeChemical) = elementalscheme(product)
detectedchemical(precursor::AbstractChemical, product::CompleteSchema) = AdductIon(precursor, product, 1) 
detectedchemical(precursor::AbstractChemical, product::AbstractScheme) = AdductIon(precursor, completescheme(precursor, product), 1) 
detectedchemical(precursor::AbstractChemical, product::StructuralChemicalScheme) = elementalscheme(precursor, product)

detectedchemical(precursor::Nothing, product::AbstractChemical) = product
detectedchemical(precursor::Nothing, product::CompleteSchemeChemical) = elementalscheme(product)
detectedchemical(precursor::Nothing, product::CompleteSchema) = throw(ArgumentError("`detectedchemical` requires precursor for scheme product."))
detectedchemical(precursor::Nothing, product::AbstractScheme) = throw(ArgumentError("`detectedchemical` requires precursor for scheme product."))
detectedchemical(precursor::Nothing, product::StructuralChemicalScheme) = throw(ArgumentError("`detectedchemical` requires precursor for scheme product."))

detectedchemical(precursor::AbstractAdductIon, product::CompleteSchemeChemical) = elementalscheme(product)
detectedchemical(precursor::T, product::CompleteSchema) where {T<:AbstractAdductIon} = throw(ArgumentError("Specific `detectedchemical(precursor::$T, scheme::CompleteSchema)` method has to be implemented."))
detectedchemical(precursor::T, product::AbstractScheme) where {T<:AbstractAdductIon} = throw(ArgumentError("Specific `detectedchemical(precursor::$T, scheme::AbstractScheme)` method has to be implemented."))
detectedchemical(precursor::AbstractAdductIon, product::StructuralChemicalScheme) = elementalscheme(precursor, product)

detectedchemical(precursor::AdductIon, product::CompleteSchemeChemical) = elementalscheme(product)
detectedchemical(precursor::AdductIon, product::CompleteSchema) = AdductIon(ioncore(precursor), adductionscheme(precursor, product), ncore(precursor))
detectedchemical(precursor::AdductIon, product::AbstractScheme) = AdductIon(ioncore(precursor), adductionscheme(precursor, product), ncore(precursor))
detectedchemical(precursor::AdductIon, product::StructuralChemicalScheme) = elementalscheme(precursor, product)

detectedchemical(precursor::Isotopomers, product::AbstractChemical) = Isotopomers(detectedchemical(chemicalparent(precursor), product), Pair{String, Int}[])
detectedchemical(precursor::Isotopomers, product::Isotopomers) = product
detectedchemical(precursor::Isotopomers, product::Groupedisotopomers) = Isotopomers(chemicalparent(product), isotopomersisotopes(product))
detectedchemical(precursor::Isotopomers, product::CompleteSchemeChemical) = detectedchemical(precursor, elementalscheme(product))
function detectedchemical(precursor::Isotopomers, product::CompleteSchema)
    chemical = detectedchemical(chemicalparent(precursor), chemicalparent(product))
    isotopes = gain_elements(isotopomersisotopes(precursor), isotopomersisotopes(product))
    Isotopomers(chemical, isotopes)
end
function detectedchemical(precursor::Isotopomers, product::AbstractScheme)
    product = completescheme(precursor, product)
    chemical = detectedchemical(chemicalparent(precursor), chemicalparent(product))
    isotopes = gain_elements(isotopomersisotopes(precursor), isotopomersisotopes(product))
    Isotopomers(chemical, isotopes)
end
detectedchemical(precursor::Isotopomers, product::StructuralChemicalScheme) = detectedchemical(precursor, elementalscheme(chemicalparent(precursor), product))

detectedchemical(precursor::Groupedisotopomers, product::AbstractChemical) = Groupedisotopomers(detectedchemical(chemicalparent(precursor), product), 0, precursor.isotope, groupedisotopomersisotopes(product), groupedisotopomersabundance(product))
detectedchemical(precursor::Groupedisotopomers, product::Isotopomers) = Groupedisotopomers(chemicalparent(product), isotopomerstate(product; isotope = precursor.isotope), precursor.isotope, groupedisotopomersisotopes(product), groupedisotopomersabundance(product))
detectedchemical(precursor::Groupedisotopomers, product::Groupedisotopomers) = product
detectedchemical(precursor::Groupedisotopomers, product::CompleteSchemeChemical) = detectedchemical(precursor, elementalscheme(product))
function detectedchemical(precursor::Groupedisotopomers, product::CompleteSchema)
    chemical = detectedchemical(chemicalparent(precursor), chemicalparent(product))
    isotopes = map((x, y) -> gain_elements(x, y), groupedisotopomersisotopes(precursor), groupedisotopomersisotopes(product))
    state = _isotopomerstate(first(isotopes), elements_mass()[precursor.isotope] - elements_mass()[elements_parents()[precursor.isotope]])
    Groupedisotopomers(chemical, state, precursor.isotope, groupedisotopomersisotopes(isotopes), groupedisotopomersabundance(product))
end
function detectedchemical(precursor::Groupedisotopomers, product::AbstractScheme)
    product = completescheme(precursor, product)
    chemical = detectedchemical(chemicalparent(precursor), chemicalparent(product))
    isotopes = map((x, y) -> gain_elements(x, y), groupedisotopomersisotopes(precursor), groupedisotopomersisotopes(product))
    state = _isotopomerstate(first(isotopes), elements_mass()[precursor.isotope] - elements_mass()[elements_parents()[precursor.isotope]])
    Groupedisotopomers(chemical, state, precursor.isotope, groupedisotopomersisotopes(isotopes), groupedisotopomersabundance(product))
end
detectedchemical(precursor::Groupedisotopomers, product::StructuralChemicalScheme) = detectedchemical(precursor, elementalscheme(chemicalparent(precursor), product))

chemicalentity(isobars::Isobars; kwargs...) = chemicalentity(first(chemicalspecies(isobars)))
chemicalentity(isotopomers::Groupedisotopomers; kwargs...) = Isotopomers(chemicalparent(isotopomers), isotopomersisotopes(isotopomers))
chemicalentity(ct::ChemicalTransition; kwargs...) = chemicalentity(first(chemicaltransition(ct)))

elementalscheme(sch::Groupedisotopomerizedschema; kwargs...) = Groupedisotopomerizedschema(elementalscheme(sch.parent; kwargs...), sch.state, sch.isotope, sch.isotopes, sch.abundance)
elementalscheme(sch::IsotopomerizedSchema; kwargs...) = IsotopomerizedSchema(elementalscheme(sch.parent; kwargs...), sch.isotopes)
elementalscheme(sch::ChemicalSchema; kwargs...) = ChemicalSchema(elementalscheme.(sch.schema; kwargs...), sch.number)
structuralscheme(sch::Groupedisotopomerizedschema; kwargs...) = Groupedisotopomerizedschema(structuralscheme(sch.parent; kwargs...), sch.state, sch.isotope, sch.isotopes, sch.abundance)
structuralscheme(sch::IsotopomerizedSchema; kwargs...) = IsotopomerizedSchema(structuralscheme(sch.parent; kwargs...), sch.isotopes)
structuralscheme(sch::ChemicalSchema; kwargs...) = ChemicalSchema(structuralscheme.(sch.schema; kwargs...), sch.number)
structuralscheme(::Nothing; kwargs...) = nothing 
structuralscheme(x::Symbol; kwargs...) = x 
elementalscheme(::Nothing; kwargs...) = nothing 
elementalscheme(x::Symbol; kwargs...) = x 

chemicalspecies(isobars::Isobars; kwargs...) = isobars.chemicals

function chemicaltransition(isobars::Isobars{<:ChemicalTransition}; kwargs...) 
    ct = chemicaltransition.(chemicalspecies(isobars))
    [Isobars(getindex.(ct, i), abundance) for (i, abundance) in enumerate(eachcol(isobars.abundance))]
end
chemicaltransition(ct::ChemicalTransition; kwargs...) = ct.transition

chemicalparent(isobars::Isobars; kwargs...) = chemicalparent(chemicalentity(isobars); kwargs...)
chemicalparent(isotopomers::Isotopomers; kwargs...) = isotopomers.parent 
chemicalparent(isotopomers::Groupedisotopomers; kwargs...) = isotopomers.parent 
chemicalparent(ct::ChemicalTransition; kwargs...) = ChemicalTransition(chemicalparent.(chemicaltransition(ct); kwargs...))

chemicalparent(sch::StructuralElementalScheme; kwargs...) = StructuralElementalScheme(structuralscheme(sch), chemicalparent(elementalscheme(sch); kwargs...))
chemicalparent(sch::ElementalScheme{T}; kwargs...) where T = ElementalScheme(T, chemicalparent(sch.chemical; kwargs...))
chemicalparent(sch::IsotopomerizedSchema; kwargs...) = sch.parent
chemicalparent(sch::ChemicalSchema; kwargs...) = ChemicalSchema(chemicalparent.(sch.schema; kwargs...), sch.number)
chemicalparent(sch::Groupedisotopomerizedschema; kwargs...) = sch.parent 

inputchemical(isobars::Isobars; kwargs...) = Isobars([inputchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
inputchemical(ct::ChemicalTransition; kwargs...) = first(chemicaltransition(ct))

outputchemical(isobars::Isobars; kwargs...) = Isobars([outputchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
outputchemical(ct::ChemicalTransition; kwargs...) = last(chemicaltransition(ct))

analyzedchemical(isobars::Isobars; kwargs...) = detectedchemical(isobars; kwargs...)
analyzedchemical(isobars::Isobars{<:ChemicalTransition}; kwargs...) = Isobars([analyzedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, begin])
analyzedchemical(ct::ChemicalTransition; kwargs...) = detectedchemical(inputchemical(ct); kwargs...)

function seriesanalyzedchemical(isobars::Isobars{<:ChemicalTransition}; kwargs...) 
    ct = chemicaltransition.(chemicalspecies(isobars))
    [analyzedchemical(Isobars(getindex.(ct, i), abundance); kwargs...) for (i, abundance) in enumerate(eachcol(isobars.abundance))]
end
function seriesanalyzedchemical(ct::ChemicalTransition; kwargs...) 
    v = AbstractChemical[]
    precursor = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedchemical(c; precursor))
        precursor = last(v)
    end
    v
end
function seriesanalyzedisotopes(ct::ChemicalTransition; kwargs...)
    v = Vector{Pair{String, Int}}[]
    precursorisotopes = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedisotopes(c; precursorisotopes))
        precursorisotopes = last(v)
    end
    v
end
function seriesanalyzedcharge(ct::ChemicalTransition; kwargs...)
    v = Int[]
    precursorcharge = nothing
    for c in chemicaltransition(ct) 
        push!(v, detectedcharge(c; precursorcharge))
        precursorcharge = last(v)
    end
    v
end

detectedchemical(isobars::Isobars; kwargs...) = Isobars([detectedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance)
detectedchemical(isobars::Isobars{<:ChemicalTransition}; kwargs...) = Isobars([detectedchemical(chemical; kwargs...) for chemical in chemicalspecies(isobars)], isobars.abundance[:, end])

detectedisotopes(sch::CompleteSchemeChemical; precursor = nothing, precursorisotopes = nothing, kwargs...) = isotopomersisotopes(sch; kwargs...)
detectedcharge(sch::CompleteSchemeChemical; precursor = nothing, precursorcharge = nothing, kwargs...) = charge(sch; kwargs...)
detectedelements(sch::CompleteSchemeChemical; precursor = nothing, precursorelements = nothing, kwargs...) = chemicalelements(sch; kwargs...)
# detectedchemical(ct::ChemicalTransition; kwargs...) = last(seriesanalyzedchemical(ct; kwargs...))
# detectedisotopes(ct::ChemicalTransition; kwargs...) = last(seriesanalyzedisotopes(ct; kwargs...))::Vector{Pair{String, Int}}
# detectedcharge(ct::ChemicalTransition; kwargs...) = last(seriesanalyzedcharge(ct; kwargs...))::Int
# detectedelements(ct::ChemicalTransition; kwargs...) = last(seriesanalyzedelements(ct; kwargs...))::Vector{Pair{String, Int}}
# detectedchemical(ct::ChemicalTransition; kwargs...) = 
#     detectedchemical(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...) 
# detectedchemical(ct::AbstractVector; kwargs...) = 
#     length(ct) > 1 ? detectedchemical(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) : first(ct)
    # detectedchemical(@view(chemicaltransition(ct)[begin:end - 1]), outputchemical(ct); kwargs...) 
# detectedchemical(ct::AbstractVector, chemical::AbstractChemical; kwargs...) = chemical
function detectedchemical(ct::ChemicalTransition; kwargs...) 
    precursor = nothing
    for c in chemicaltransition(ct) 
        precursor = detectedchemical(precursor, c)
    end
    precursor
end
detectedisotopes(ct::ChemicalTransition; kwargs...) = 
    detectedisotopes(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...)::Vector{Pair{String, Int}}
detectedisotopes(ct::AbstractVector; kwargs...) = 
    (length(ct) > 1 ? detectedisotopes(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) : isotopomersisotopes(first(ct); kwargs...))::Vector{Pair{String, Int}}
detectedcharge(ct::ChemicalTransition; kwargs...) = 
    detectedcharge(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...)::Int
detectedcharge(ct::AbstractVector; kwargs...) = 
    (length(ct) > 1 ? detectedcharge(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) : charge(first(ct); kwargs...))::Int
# detectedelements(ct::ChemicalTransition; kwargs...) = 
#     detectedelements(outputchemical(ct); precursor = @view(chemicaltransition(ct)[begin:end - 1]), kwargs...)::Vector{Pair{String, Int}}
# detectedelements(ct::AbstractVector; kwargs...) = 
#     (length(ct) > 1 ? detectedelements(last(ct); precursor = @view(ct[begin:end - 1]), kwargs...) : chemicalelements(first(ct); kwargs...))::Vector{Pair{String, Int}}

# detectedelements(ct::ChemicalTransition; kwargs...) = chemicalelements(detectedchemical(ct); kwargs...)

# GenericChemical property search
# structural -> complete, elemental -> complete
structure_search(chemical, precursor_schema, product_schema::CompleteSchema) = product_schema
structure_search_elemental(chemical, precursor_schema, product_schema::IsotopomerizedSchema) = IsotopomerizedSchema(structure_search_elemental(chemical, precursor_schema, product_schema.parent), product_schema.isotopes)
structure_search(chemical, precursor_schema, product_schema::IsotopomerizedSchema) = IsotopomerizedSchema(structure_search(chemical, precursor_schema, product_schema.parent), product_schema.isotopes)
structure_search_elemental(chemical, precursor_schema, product_schema::ChemicalSchema) = ChemicalSchema([structure_search_elemental(chemical, precursor_schema, k) for k in product_schema.schema], product_schema.number)
structure_search(chemical, precursor_schema, product_schema::ChemicalSchema) = ChemicalSchema([structure_search(chemical, precursor_schema, k) for k in product_schema.schema], product_schema.number)
function structure_search(chemical, precursor_schema, product::GenericChemical) 
    product_schema = chemicalscheme_search(product, nothing) 
    isnothing(product_schema) && return StructuralElementalScheme(RandomProductScheme(), product)
    sch = _structure_search(chemical, precursor_schema, product_schema)
    isnothing(sch) ? StructuralElementalScheme(RandomProductScheme(), product) : StructuralElementalScheme(structuralscheme(product_schema), sch)
end
function structure_search(chemical, precursor_schema, product::AdductIon{<:GenericChemical}) 
    product_schema = chemicalscheme_search(ioncore(product), ionadduct(product)) 
    isnothing(product_schema) && return StructuralElementalScheme(RandomProductScheme(), product)
    sch = _structure_search(chemical, precursor_schema, product_schema)
    isnothing(sch) ? StructuralElementalScheme(RandomProductScheme(), product) : StructuralElementalScheme(structuralscheme(product_schema), sch)
end
function structure_search_elemental(chemical, precursor_schema, product::GenericChemical) 
    product_schema = chemicalscheme_search(product, nothing) 
    isnothing(product_schema) && return product
    sch = _structure_search(chemical, precursor_schema, product_schema)
    isnothing(sch) ? product : sch
end
function structure_search_elemental(chemical, precursor_schema, product::AdductIon{<:GenericChemical}) 
    product_schema = chemicalscheme_search(ioncore(product), ionadduct(product)) 
    isnothing(product_schema) && return product
    sch = _structure_search(chemical, precursor_schema, product_schema)
    isnothing(sch) ? product : sch
end

function structure_search(chemical, precursor_schema, product_schema::AbstractElementalScheme) 
    StructuralElementalScheme(product_schema, structure_search_elemental(chemical, precursor_schema, product_schema))
end

function structure_search_elemental(chemical, precursor_schema, product_schema::AbstractElementalScheme) 
    sch = _structure_search(chemical, precursor_schema, product_schema)
    isnothing(sch) ? copy(product_schema) : sch
end

function _structure_search(chemical, precursor_schema, product_schema)
    schema = getchemicalproperty(chemical, :structure, nothing)
    isnothing(schema) && return nothing
    i = findfirst(x -> first(x) == structuralscheme(precursor_schema), schema)
    isnothing(i) && return nothing
    scheme = last(schema[i])
    i = findfirst(x -> first(x) == structuralscheme(product_schema), scheme)
    isnothing(i) && return nothing
    last(scheme[i])
end

function chemicalscheme_search(chemical, scheme)
    schema = getchemicalproperty(chemical, :chemicalscheme, nothing)
    isnothing(schema) && return nothing
    i = findfirst(x -> first(x) == structuralscheme(scheme), schema)
    isnothing(i) && return nothing
    last(schema[i])
end

# scheme -> complete
schema_search(chemical, precursor_schema::Nothing, product_schema::CompleteSchema) = product_schema
schema_search(chemical, precursor_schema::Nothing, product_schema) = completescheme(chemical, product_schema)
function schema_search(chemical, precursor_schema, product_schema)
    schema = getchemicalproperty(chemical, :schema, nothing)
    isnothing(schema) && return nothing
    i = findfirst(x -> first(x) == structuralscheme(precursor_schema), schema)
    # isnothing(i) && return StructuralElementalScheme(ChemicalSchema(structuralscheme(precursor_schema), structuralscheme(product_schema)), ChemicalSchema(elementalscheme(precursor_schema), elementalscheme(product_schema)))
    isnothing(i) && return nothing
    scheme = last(schema[i])
    i = findfirst(x -> first(x) == structuralscheme(product_schema), scheme)
    # isnothing(i) && return StructuralElementalScheme(ChemicalSchema(structuralscheme(precursor_schema), structuralscheme(product_schema)), ChemicalSchema(elementalscheme(precursor_schema), elementalscheme(product_schema)))
    isnothing(i) && return nothing
    last(scheme[i])
end

