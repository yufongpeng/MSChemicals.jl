module MassSpecChemicals

using Combinatorics, TypedTables, MLStyle, Statistics, StatsBase, Dictionaries, Intervals, SplitApplyCombine, Plots
using UnitfulMoles: parse_compound
using SentinelArrays: ChainedVector
import Base: show, length, +, -, *, /, isless, isequal, in, union, intersect, iterate, Broadcast.broadcastable, ==, hash, copy, axes, isempty, eltype

export 
    # Types
    AbstractChemical, Chemical, FormulaChemical, ChemicalTransition, ChemicalSeries, Isobars, Isotopomers, AbstractChemicalWrapper, 
    AbstractAdductIon, AdductIon,
    AbstractScheme, AbstractElementalScheme, AbstractStructuralScheme, StructuralChemicalScheme, AbstractCompleteScheme, 
    ElementalScheme, ChemicalGain, ChemicalLoss, 
    ChemicalSchema, IsotopomerizedSchema, StructuralElementalScheme, 
    CompleteSchema, StructuralSchema, ElementalSchema, 

    # Default chemical for scheme 
    Electron, 
    Proton, 
    Water, 
    Ammonia, 
    Ammonium, 
    Sodium, 
    Potassium, 
    Lithium, 
    Silver, 
    Acetate, 
    Formate, 
    AceticAcid, 
    FormicAcid, 
    MethylAcetate, 
    MethylFormate, 
    Fluoride, 
    Chloride, 
    Methenium, 
    
    # Setter 
    set_scheme!, set_schabbr!, set_element!, 

    # Elements 
    parent_element, major_isotope, minor_isotope, iselement, isisotope, 

    # Parser
    parse_chemical, parse_adduct, 
    ChemicalParser, FormulaChemicalParser, AdductParser, ChemicalExpressionParser, ChemicalTransitionParser, ChemicalEntityParser, ChemicalSchemeParser, 

    # Attributes
    chemicalformula, chemicalelements, chemicalname, chemicalabbr, chemicalsmiles, ioncore, ionadduct, ncore, charge, ncharge, retentiontime,
    mmi, molarmass, mz, 
    chemicalentity, chemicalspecies, chemicalpair, 
    chemicalparent, isotopomersisotopes, groupedisotopomersisotopes, groupedisotopomersabundance,
    analyzedchemical, detectedchemical, detectedelements, detectedcharge, detectedisotopes,
    inputchemical, outputchemical,
    seriesanalyzedchemical, seriesanalyzedisotopes, seriesanalyzedcharge, seriesanalyzedelements, 
    msstage, chemicaltransition, isotopomerstate, 
    getchemicalproperty, 
    elementalscheme, structuralscheme, 

    # Scheme 
    isgainscheme, islossscheme, completescheme, adductionscheme, isotopomerize, ionize, 

    # Derived attributes
    isotopicabundance, Isotopologues, TandemIsotopologues, 
    CoelutingIsobars, isobar_table, group_isotopologues, 
    
    # Spectrum 
    Ionization, Spectrum, MSScan, AllIons, Isolation, SelectedIonMonitor, Fragmentation, peak_table, 
    MSAnalyzer, Quadrupole, QuadrupoleIonTrap, QIT, LinearIonTrap, LIT, TimeOfFlight, TOF, Orbitrap, FourierTransformIonCyclotronResonance, FTICR, 
    resolving_power, 
    
    # Visualization
    plot_spectrum, plot_spectrum!, 
    plot_resolving_power, plot_resolving_power!,
    plot_window, plot_window!,

    # Utils
    match_chemical, 
    ischemicalequal, isadductequal, 
    acrit, rcrit, crit, @ri_str, makecrit_delta, makecrit_value, 
    value_error, relative_error, percentage_error, ppm_error, relative_error_mean, percentage_error_mean, ppm_error_mean


abstract type AbstractChemicalsSchema end
"""
    AbstractChemical <: AbstractChemicalsSchema

Abstract type for chemicals. 
    
The attribute function `chemicalname -> String` (unique chemical name) is required for a concrete type of `AbstractChemical`. 
It defaults to property `name`.

At least one of the following attributes are required. They are interchangable.
* `chemicalformula -> String`: chemical formula. It defaults to property `:formula`.
* `chemicalelements -> Vector{Pair{String, Int}}`: chemical elements. It defaults to property `:elements`.

The following attributes are optional, but generic functions are defined.
* `chemicalabbr -> String`: abbreviation. It defaults to property `:abbreviation` and `chemicalname`. 
* `chemicalsmiles -> String`: SMILES. It defaults to property `:SMILES` and `""`.
* `charge -> Int`: charge state; positive for cation, negative for anion. It defaults to property `:charge` annd `0`.
* `ncharge -> Int`: number of charges.
* `retentiontime -> AbstractFloat`: retention time. It defaults to property `:retentiontime` and `NaN`.
* `mmi -> AbstractFloat`: monoisotopic mass.
* `molarmass -> AbstractFloat`: molar mass.
* `mz -> AbstractFloat`: mass to charge ratio (m/z).
* `chemicalparent -> AbstractChemical`: parent chemical without delocalized isotopes replacement.
* `isotopomersisotopes -> Vector{Pair{String, Int}}`: delocalized isotopes replacement of isotopomers.
* `isotopomerstate -> Int`: isotopomers state, i.e. equivalent number of isotope.
* `groupedisotopomersisotopes -> Vector{Vector{Pair{String, Int}}}`: delocalized isotopes replacements of each isotopomers.
* `groupedisotopomersabundance -> AbstractFloat`: abundance of each isotopomers.
* `chemicalentity -> AbstractChemical`: a single chemical entity representing the chemical.
* `chemicalspecies -> Vector{<: AbstractChemical}`: multiple chemical entities having shared properties. 
* `chemicaltransition -> Vector{<: AbstractChemical}`: chemical entities analyzed in each stage of instrumental analysis.  
* `inputchemical -> AbstractChemical`: a single chemical entity that is the input at the very beginning of instrumental analysis. 
* `outputchemical -> AbstractChemical`: a single chemical entity that is the output at the very end of instrumental analysis. 
* `analyzedchemical -> AbstractChemical`: a single chemical entity directly detected at the very beginning of instrumental analysis. 
* `detectedchemical -> AbstractChemical`: a single chemical entity directly detected at the very end of instrumental analysis. 
* `seriesanalyzedchemical -> Vector{<: AbstractChemical}`: chemical entities directly analyzed in each stage of instrumental analysis.
* `detectedcharge -> Int`: charge of detected chemical.
* `detectedelements -> Vector{Pair{String, Int}}`: elements of detected chemical.
* `detectedisotopes -> Vector{Pair{String, Int}}`: isotopes replacement of detected chemical.
* `seriesanalyzedcharge -> Vector{Int}`: charge of sereially analyzed chemicals.
* `seriesanalyzedelements -> Vector{Vector{Pair{String, Int}}}`: elements of sereially analyzed chemicals.
* `seriesanalyzedisotopes -> Vector{Vector{Pair{String, Int}}}`: isotopes replacements of sereially analyzed chemicals.
* `msstage -> Int`: number of stages of MS the chemical has been through.

Specific Methods for the attributes are defined for other intrinsic chemical type on different chemical level.
* Entity Level: attribute of the corresponding chemical entity.
* Species Level: attribute of the corresponding chemical species.
* Transition Level: attribute of the corresponding chemical transition.
"""
abstract type AbstractChemical <: AbstractChemicalsSchema end

"""
    AbstractScheme <: AbstractChemicalsSchema

Abstract type for all kinds of chemical schema.

The following atributes are implemented.
* `elementalscheme -> AbstractScheme`.
* `structuralalscheme -> AbstractScheme`.
* `chemicalname -> String`.
* `chemicalformula -> String`.
* `chemicalelements -> Vector{Pair{String, Int}}`.
* `chemicalabbr -> String`: abbreviation.
* `charge -> Int`: charge state; positive for cation, negative for anion.
* `ncharge -> Int`: number of charges.
* `mmi -> AbstractFloat`: monoisotopic mass.
* `molarmass -> AbstractFloat`: molar mass.
* `chemicalparent -> AbstractChemical`: parent scheme without delocalized isotopes replacement.
* `isotopomersisotopes -> Vector{Pair{String, Int}}`: delocalized isotopes replacement of isotopomers.
* `isotopomerstate -> Int`: isotopomers state, i.e. equivalent number of isotope.
* `groupedisotopomersisotopes -> Vector{Vector{Pair{String, Int}}}`: delocalized isotopes replacements of each isotopomers.
* `groupedisotopomersabundance -> AbstractFloat`: abundance of each isotopomers.

Specific Methods for the attributes are defined for other intrinsic scheme type on different scheme level.
* Entity Level: attribute of the elemental scheme.
* Species Level: attribute of the scheme itself.
* Transition Level: only apply to schema in `ChemicalTransition`; `Species Level` for each scheme.
"""
abstract type AbstractScheme <: AbstractChemicalsSchema end
# mt, ccs 
include(joinpath("type", "chemical.jl"))
include(joinpath("type", "scheme.jl"))
include(joinpath("type", "parser.jl"))
include(joinpath("type", "utils.jl"))
include(joinpath("type", "spectrum.jl"))
include("interface.jl")
include("attr.jl")
# include("attr_scheme.jl")
include("attr_measure.jl")
include("chemical_scheme.jl")
include("elements.jl")
include("scheme.jl")
include("formula_elements.jl")
include("measure.jl")
include("isotopologues_utils.jl")
include("isotopologues_ms1.jl")
include("isotopologues_ms2.jl")
include("isotopologues_msn.jl")
include("isotopologues.jl")
include("spectrum.jl")
include("isobars.jl")
include("plot.jl")
include("input.jl")
include("output.jl")
include("utils.jl")

end