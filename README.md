# MassSpecChemicals

[![Build Status](https://github.com/yufongpeng/MassSpecChemicals.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/yufongpeng/MassSpecChemicals.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/yufongpeng/MassSpecChemicals.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yufongpeng/MassSpecChemicals.jl)

A Julia package for representing molecules and ions formed in mass spectrometers (MS).

All chemicals are instances of the abstract type `AbstractChemical`.
Charged chemicals formed in MS with a specific adduct or neutral loss (adduct ions) are instances of the abstract type `AbstractAdductIon`.

# Built-in chemical types

1. `Chemical`: unstructured chemicals storing name, elements, and other attributes

    ```julia
    Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}; property...)

    Chemical(name::AbstractString, formula::String; property...)

    Chemical(name::AbstractString, elements::Vector{Pair{String, Int}}, property::Vector{Pair{Symbol, Any}})
    ```

2. `FormulaChemical`: unstructured chemicals using the formula as the name

    ```julia
    FormulaChemical(elements::Vector{Pair{String, Int}}; property...)

    FormulaChemical(formula::AbstractString; property...)

    FormulaChemical(elements::Vector{Pair{String, Int}}, property::Vector{Pair{Symbol, Any}})
    ```

3. `AdductIon`: charged chemicals with a specific adduct or neutral loss

    ```julia
    AdductIon(core::AbstractChemical, adduct::AbstractScheme, ncore::Int)

    AdductIon(core::AbstractChemical, adduct::AbstractString)
    ```

4. `ChemicalTransition`: MS/MS transitions

    ```julia
    ChemicalTransition(transition::Vector)

    ChemicalTransition(precursor::AbstractChemical, products...)
    ```

5. `Isobars`: multiple chemicals with similar m/z

    ```julia
    Isobars(chemical::Vector{<: AbstractChemical}, abundance::VecOrMat)
    ```

6. `Isotopomers`: multiple chemicals differing by isotopic replacement location

    ```julia
    Isotopomers(parent::AbstractChemical, isotopes::ElementsVector)

    Isotopomers(parent::AbstractChemical, fullformula::AbstractString)

    Isotopomers(parent::AbstractChemical, fullelements::Dictionary)

    Isotopomers(parent::AbstractChemical, fullelements::Vector{Pair{String, Int}})
    ```

7. `Groupedisotopomers`: isotopomers grouped by isotopomer state

    ```julia
    Groupedisotopomers(parent::AbstractChemical, state::Int, isotope::String, isotopes::Vector{Vector{Pair{String, Int}}}, abundance::Vector)
    ```

Users can parse chemical expressions and pairs using `parse_chemical`.

```julia
parse_chemical("[H3O]+") # H3O with one positive charge
parse_chemical("[H2O+H]+") # protonated H2O
parse_chemical("[C6H12O6+H]+ -> [H2O+H]+") # protonated C6H12O6 fragmented to protonated H2O
parse_chemical("[C6H12O6+H]+" => "-H2O") # protonated C6H12O6 and neutral loss H2O
parse_chemical("[C6H12O6+2H]2+" => "-[H2O+H]+") # diprotonated C6H12O6 and loss of protonated H2O
parse_chemical(ChemicalTransitionParser(ChemicalExpressionParser(; charge = 2, loss = 1)), "C6H14O6" => "-H3O") # C6H14O6 (charge = 2) and loss H3O (charge = 1)
parse_chemical("[C6H12O6+2H]2+" => ChemicalLoss(Water())) # mixed parser with another chemical type
```

# Elements
## Parent elements and major isotopes

|Symbol|Major isotopes|Atomic number|Mass number|
|--------|-------------|-----------|------------------|
|C|[12C]|6|12|
|H|[1H]|1|1|
|O|[16O]|8|16|
|N|[14N]|7|14|
|P|[31P]|15|31|
|S|[32S]|16|32|
|Li|[7Li]|3|7|
|Na|[23Na]|11|23|
|K|[39K]|19|39|
|F|[19F]|9|19|
|Cl|[35Cl]|17|35|
|Ag|[108Ag]|47|108|
|Se|[80Se]|34|80|

## Minor isotopes

|Minor isotopes|Atomic number|Mass number|Alternative symbol|
|--------|-------------|-----------|------------------|
|[13C]|6|13||
|D|1|2|[2H]|
|[17O]|8|17||
|[18O]|8|18||
|[15N]|7|15||
|[33S]|16|33||
|[34S]|16|34||
|[36S]|16|36||
|[6Li]|3|6||
|[40K]|19|40||
|[41K]|19|41||
|[37Cl]|17|37||
|[109Ag]|47|109||
|[74Se]|34|74||
|[76Se]|34|76||
|[77Se]|34|77||
|[78Se]|34|78||
|[82Se]|34|82||

By default, parent elements are considered major isotopes and may be replaced by minor isotopes. For example:
* `CO2` has carbon-12 and two oxygen-16 atoms, but minor isotope replacements are possible.
* `[13C][16O]O` has carbon-13, and two oxygen-16 atoms. Carbon-13 atom and one oxygen-16 atom cannot be replaced by other isotopes.

One exception is that for parent chemical of `Isotopomers`, parent elements are major isotopes and the number of replacements is restricted by the `isotopes` field.

Users can use:

```julia
set_elements!(element, mass, abundance; minor_name = nothing)
```

to add new elements with mass and natural abundances for all isotopes.
Customized minor element names (`minor_name`) are optional.

# AbstractScheme

Any chemical gain, loss, and fragmentation scheme is an instance of `AbstractScheme`. This type has three abstract subtypes:

1. `AbstractElementalScheme`: a scheme containing elemental information, including isotopic replacement.
2. `AbstractStructuralScheme`: a scheme containing only structural information. This is useful for defining rule-based fragmentation.
3. `AbstractCompleteScheme`: a scheme containing both elemental and structural information. It is the final scheme stored in `AdductIon`.

In addition to single scheme, multiple schema are wrapped in `ChemicalSchema`.

Predefined chemicals used in adducts:

|Chemical|Abbreviation|
|-------------|-----------------|
|`Electron`|`"e"`|
|`Proton`|`"H"`|
|`Water`|`"H2O"`|
|`Ammonia`|`"NH3"`|
|`Ammonium`|`"[NH4]+"`|
|`Sodium`|`"[Na]+"`|
|`Potassium`|`"[K]+"`|
|`Silver`|`"[Ag]+"`|
|`OAc`|`"[OAc]-"`|
|`OFo`|`"[OFo]-"`|
|`Fluoride`|`"[F]-"`|
|`Chloride`|`"[Cl]-"`|
|`Methinium`|`"[Me]+"`|
|`AceticAcid`|`"HOAc"`|
|`FormicAcid`|`"HOFo"`|
|`Methylacetate`|`"MeOAc"`|
|`Methylformate`|`"MeOFo"`|

# API
## Attributes of `AbstractChemical`

Attributes are interfaces for directly accessing properties and fields through `getchemicalproperty`, or for deriving values from other attributes.

|Attribute|Return type|Description|
|----|----|-----------|
|`chemicalname`|`String`|unique chemical name|
|`chemicalformula`|`String`|chemical formula|
|`chemicalelements`|`Vector{Pair{String, Int}}`|chemical elements|
|`chemicalabbr`|`String`|common abbreviation; defaults to `chemicalname`|
|`chemicalsmiles`|`String`|SMILES; defaults to `""`|
|`charge`|`Int`|net charge (positive or negative); defaults to 0|
|`ncharge`|`Int`|number of charges|
|`retentiontime`|`Float64`|retention time; defaults to `NaN`|
|`chemicalparent`|`AbstractChemical`|parent chemical without delocalized isotope replacements|
|`isotopomersisotopes`|`Vector{Pair{String, Int}}`|delocalized isotope replacements of isotopomers|
|`isotopomerstate`|`Int`|isotopomer state, i.e. equivalent number of isotopes|
|`groupedisotopomersisotopes`|`Vector{Vector{Pair{String, Int}}}`|delocalized isotope replacements of each isotopomer in group|
|`groupedisotopomersabundance`|`AbstractFloat`|abundance of each isotopomer in group|
|`chemicalentity`|`AbstractChemical`|a single chemical entity representing the chemical|
|`chemicalspecies`|`Vector{<: AbstractChemical}`|multiple chemical entities with shared properties|
|`chemicaltransitions`|`Vector{<: AbstractChemical}`|chemical entities analyzed in each stage of instrumental analysis|
|`inputchemical`|`AbstractChemical`|the chemical entity that is the input at the beginning of analysis|
|`outputchemical`|`AbstractChemical`|the chemical entity that is the output at the end of analysis|
|`analyzedchemical`|`AbstractChemical`|the chemical entity directly analyzed at the beginning of analysis|
|`detectedchemical`|`AbstractChemical`|the chemical entity directly detected at the end of analysis|
|`detectedisotopes`|`Vector{Pair{String, Int}}`|delocalized isotope replacements of the detected chemical|
|`detectedcharge`|`Int`|charge state of the detected chemical|
|`detectedelements`|`Vector{Pair{String, Int}}`|elements of the detected chemical|
|`seriesanalyzedchemical`|`Vector{<: AbstractChemical}`|chemical entities directly analyzed in each stage of instrumental analysis|
|`seriesanalyzedisotopes`|`Vector{Vector{Pair{String, Int}}}`|delocalized isotope replacements of serially analyzed chemicals|
|`seriesanalyzedcharge`|`Vector{Int}`|charge states of serially analyzed chemicals|
|`seriesanalyzedelements`|`Vector{Vector{Pair{String, Int}}}`|elements of serially analyzed chemicals|
|`msstage`|`Int`|number of MS stages the chemical has been through|
|`mmi`|`Float64`|monoisotopic mass|
|`molarmass`|`Float64`|molar mass|
|`mz`|`Float64`|m/z, mass-to-charge ratio|

Specific methods for attributes are defined for each intrinsic chemical type at different chemical levels:
* Entity Level: attribute of the corresponding chemical entity
* Species Level: attribute of the corresponding chemical species
* Transition Level: attribute of the corresponding chemical transitions

See documentation for each attribute to determine the exact level.

### Type-specific attributes

Users can define new attributes or overload existing attribute functions for any chemical type.

```julia
abstract type Lipid <: AbstractChemical end
struct Fattyacid <: Lipid
    ncarbon::Int
    ndoublebond::Int
end
struct Acylglycerol <: Lipid
    carbonchains::Vector{Fattyacid}
end

# Different attribute methods for different chemical types
ncarbon(chemical::Fattyacid) = chemical.ncarbon
ncarbon(chemical::Acylglycerol) = sum(ncarbon, chemical.carbonchains) + 3

fa1 = Fattyacid(18, 0)
fa2 = Fattyacid(18, 1)
dg = Acylglycerol([fa1, fa2])

ncarbon(fa1) == 18
ncarbon(dg) == 39
```

### Instance-specific attributes

For attributes defined in only some instances, users can use `getchemicalproperty` to handle missing attributes.

```julia
# Access the :carbonchains property through getchemicalproperty, defaulting to [chemical]
carbonchains(chemical::Lipid) = getchemicalproperty(chemical, :carbonchains, [chemical])
carbonchains(dg) == [fa1, fa2]
carbonchains(fa1) == [fa1]
```

The type `Chemical` stores any non-default attributes in the field `property`. Users can create objects with these attributes as keyword arguments, or add attributes by directly pushing the `attr_name => attr_value` pair to `chemical.property`, then defining an attribute function to access it.

```julia
chemical = Chemical("18:0 PC-d9", "	C44H79NO8PD9"; lipidclass = "PC")

lipidclass(chemical::Chemical) = getchemicalproperty(chemical, :lipidclass) # new attribute
lipidclass(chemical) == "PC"

push!(chemical.property, :retentiontime => 10)
retentiontime(chemical) == 10 # Defined as accessing :retentiontime through getchemicalproperty
```

### Additional Atttributes of `AbstractAdductIon`
|Attribute|Return type|Description|
|----|----|-----------|
|`ioncore`|`AbstractChemical`|core chemical|
|`ionadduct`|`AbstractAdduct`|adduct originated from ionization|
|`ncore`|`Int`|number of core chemical|

When isotopes are involved in addut ion formation for an object `adduct_ion` which `chemical = ioncore(adduct_ion)::ChemicalType` and `adduct = ionadduct(adduct_ion)::Existing_Scheme`, there are two solutions.
1. If `ChemicalType` is a customized chemical type, define type-specific `completescheme`
    ```julia
    elementalscheme(adduct_ion::ChemicalType, ::Affected_Scheme) # Ionization
    elementalscheme(adduct_ion::AdductIon{ChemicalType, Existing_Scheme}, ::Affected_Scheme) # Fragmentation (Neutral Loss)
    ```
    
    `elementalscheme` returns elemental scheme containing correct elements. For instance, [M-Me]- of Deuterium-labeled phosphatidylcholine (as type `DLPC` for instance) may turn out to be [M-CD3]- (`ElementalScheme{false, DLMe}`) rather than [M-CH3]- (`ElementalScheme{false, Me}`) if Deuteriums are labeled on the methyl group of choline. In this case, extend `elementalscheme(::DLPC, ::ElementalScheme{false, Me})` such that
    ```julia 
    struct PC <: AbstractChemical 
        ... 
    end
    struct DLPC <: AbstractChemical 
        ... 
    end
    struct Me <: AbstractChemical 
        ... 
    end
    struct DLMe <: AbstractChemical 
        ... 
    end
    # dlpc: [M-Me]- of Deuterium-labeled phosphatidylcholine on the methyl group of choline
    # dlpc: [M-Me]- of Deuterium-labeled phosphatidylcholine on other part
    # pc: [M-Me]- of natural phosphatidylcholine 
    # loss_me: [M-CH3]-, i.e. ElementalScheme{false, Me}
    # loss_cd3: [M-CD3]-, i.e. ElementalScheme{false, DLMe}
    elementalscheme(dlmcpc, loss_me) == loss_cd3
    elementalscheme(dlpc, loss_me) == loss_me
    elementalscheme(pc, loss_me) == loss_me
    ```
    For more details, see example in file `test/objects/customized.jl`.
2. If `ChemicalType` is `Chemical`, define an attribute `:structure` for the `chemical`. The attribute should be ionadduct-(scheme-scheme pairs) pairs. `structure_search` finds this attribute, and extracts the value of key `adduct`. 
    ```julia
    # loss_me: [M-CH3]-, i.e. ElementalScheme(false, Chemical("Me", "CH3")}
    # loss_cd3: [M-CD3]-, i.e. ElementalScheme(false, Chemical("Me[D3]", "CD3"))
    chemical = Chemical("18:0 PC-d9", "C44H79NO8PD9")
    loss_me = ElementalScheme(false, Chemical("Me", "CH3"))
    loss_cd3 = ElementalScheme(false, Chemical("Me[D3]", "CD3"))
    push!(chemical.property, :structure => [nothing => [loss_me => loss_cd3]])
    elementalscheme(chemical, ChemicalGain(Proton())) == ChemicalGain(Proton()) # No key Protonation()
    elementalscheme(chemical, loss_me) == loss_cd3
    ```
    For more details, see example in file `test/objects/generic.jl`.
## Isotopic abundance and Isotopologues
There are three related functions
* `isotopicabundance`
    
    This function computes isotopic abundance of the input elements composition (`Vector{Pair{Int}}` or `Dictionary`), formula, and chemical (converted to elements by `chemicalelements`). Parent elements are viewed as major isotopes, and isotopic abundances of all elements are considered in computation. To compute isotopic abundance of chemicals with all isotopes labeled intentionally and not following natural distribution, set keyword argument `ignore_isotopes` true, and only parent elements are considered.

* `Isotopologues`

    This function computes isotopologues of formula, single chemical, MS/MS precursor-product pairs (formula pairs or `ChemimcalTransition`) or multiple chemicals. Only isotopic abundance of parent elements are considered, and isotopes are viewed as intentionally labeled elements. Isotopologues can be filtered by abundance threshold. MS/MS product can be neutral loss/gain or chemical loss/gain. 

* `TandemIsotopologues`

    This function is similar to `Isotopologues`; it computes isotopologues of given precursor(s) and additionally computes the abundance of fragments with given fragmentation patterns. The key difference is that this function is recursive and abundance is calculated from the beginning. It generally performs slightly slower for multiple MS stages and abundance is normalized in the first stage and filtered in all stages.

Isotopologues table can be aggregated using `group_isotopologues`.
## Mass Spectrometer 
There are five functions to simulate ions in mass spectrometer. 
1. `Isolation`: isolating target ion(s) with specific m/z values and resolution to enter the next MS stage. 
2. `AllIons`: allow all Ions within m/z range entering the next MS stage. 
3. `Fragmentation`: create a table of fragments with given fragmentation patterns. 
4. `MSScan`: perform MS scan using given mass analyzer. This function creates `Spectrum` objects, which can be visualized by `plot_spectrum` and `plot_spectrum!`. Peak lists can be extracted with function `peak_table`. 
5. `SelectedIonMonitor`: perform selected ion monitoring for target transition(s). Peak lists can be extracted with function `peak_table`. 

The following common mass analyzer are defined.
1. `Quadrupole`
2. `QuadrupoleIon` or `QIT`
3. `LinearIonTrap` or `LIT`
4. `TimeOfFlight` or `TOF`
5. `Orbitrap`
6. `FourierTransformIonCyclotronResonance` or `FTICR`

Default settings related to resolution, and isolation window are also defined for each analyzer. To create generic mass analyzer, use `MSAnalyzer`.
## Co-eluting isobars
The function `CoelutingIsobars` creates an object `CoelutingIsobars` with a vector of elution function-criteria pair, a vector of ms analyzer-criteria pair, and a target chemical table.

This object can be further aggregated using `isobar_table`.

## Other Functions
* `parent_element`: return parent element.
* `major_isotope`: return major isotope.
* `minor_isotope`: return minor isotope.
* `iselement`: whether it is an element.
* `isisotope`: whether it is an isotope (including element).
* `ischemicalequal`: whether two chemicals are chemically equivalent.
* `ischemicalequaltransform`: return an object for comparison with other chemicals by `istransformedadduct`.
* `istransformedchemicalequal`: whether two chemicals are chemically equivalent after applying `istransformedchemicalequal`.
* `isadductequal`: whether two adducts chemically equivalent.
* `isadductequaltransform`: return an object for comparison with other adducts by `istransformedadduct`.
* `istransformedadduct`: whether two adducts are chemically equivalent after applying `isadductequaltransform`.
* `match_chemical`: match detected chemicals with reference library.
* `acrit`: create absolute criterion.
* `rcrit`: create relative criterion.
* `crit`: create both absolute and relative criterion.
* `@ri_str`: real number interval.
* `plot_resolving_power` and `plot_resolving_power!`: plot the function of m/z to resolving_power.
* `plot_window` and `plot_window!`: plot the window function.