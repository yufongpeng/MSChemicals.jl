const SCHEME_ABBR = Dict{String, AbstractChemical}(
    "H"         => Proton(), 
    "H2O"       => Water(), 
    "NH3"       => Ammonia(), 
    "NH4"       => Ammonium(), 
    "Na"        => Sodium(), 
    "K"         => Potassium(), 
    "Li"        => Lithium(), 
    "Ag"        => Silver(), 
    "CH3COO"    => Acetate(), 
    "OAc"       => Acetate(), 
    "AcO"       => Acetate(), 
    "HCOO"      => Formate(), 
    "OFo"       => Formate(), 
    "FoO"       => Formate(), 
    "F"         => Fluoride(),
    "Cl"        => Chloride(),
    "Me"        => Methenium(),
    "HOAc"      => AceticAcid(),
    "AcOH"      => AceticAcid(),
    "AcA"       => AceticAcid(),
    "CH3COOH"   => AceticAcid(),
    "HOFo"      => FormicAcid(),
    "FoOH"      => FormicAcid(),
    "FoA"       => FormicAcid(),
    "FA"        => FormicAcid(),
    "HCOOH"     => FormicAcid(),
    "MeOAc"     => MethylAcetate(),
    "CH3OAc"    => MethylAcetate(),
    "AcOMe"     => MethylAcetate(),
    "AcOCH3"    => MethylAcetate(),
    "CH3COOCH3" => MethylAcetate(),
    "MeOFo"     => MethylFormate(),
    "CH3OFo"    => MethylFormate(),
    "FoOMe"     => MethylFormate(),
    "FoOCH3"    => MethylFormate(),
    "HCOOCH3"   => MethylFormate()
)

const SCHEME_NAME = Dict{String, AbstractScheme}()

scheme_doc = """
    scheme_name()
    scheme_abbr()

Access constants related to scheme. 
* `scheme_name`: scheme type for an scheme expression.
* `scheme_abbr`: conversion of abbreviation to chemical for scheme.
"""
@doc scheme_doc
scheme_name() = SCHEME_NAME

@doc scheme_doc
scheme_abbr() = SCHEME_ABBR

"""
    set_schabbr!(abbr::AbstractString, chemical::AbstractChemical)

Set `abbr` to be an abbreviation of scheme chemical `fm` for default [`AdductParser`](@ref). 

# Examples
```julia
julia> serine = Chemical("Serine", "C3H5NO2"; abbreviation = "Ser") 
Serine

julia> set_schabbr!("Ser", serine)
Dict{String, AbstractChemical} with 37 entries:
  "HCOOH"     => Formic Acid
  "Me"        => Methenium
  ⋮           => ⋮

julia> parse_adduct("[M-H-Ser]-")
(adduct = Loss_Proton|Loss_Serine, ncore = 1)
```
"""
set_schabbr!(abbr::AbstractString, chemical::AbstractChemical) = push!(scheme_abbr(), abbr => chemical)

for (abbr, fm) in [
        ("FA", FormicAcid()), 
        ("FoOH", FormicAcid()), 
        ("FoA", FormicAcid()), 
        ("HOFo", FormicAcid()), 
        ("OFo", Formate()), 
        ("FoO", Formate()), 
        ("AcA", AceticAcid()), 
        ("AcOH", AceticAcid()), 
        ("HOAc", AceticAcid()), 
        ("OAc", Acetate()), 
        ("AcO", Acetate()), 
        ("Me", Methenium())]
    set_schabbr!(abbr, fm)
end

"""
    set_scheme!(nm::AbstractString, scheme::AbstractScheme)

Set `nm` to be `scheme` for default `AdductParser`.

For customized scheme types, this function is required to make `nm` parsed into `scheme` by [`parse_adduct`](@ref).

# Examples
```julia
julia> struct NewScheme <: AbstractCompleteScheme end 

julia> set_scheme!("[+NSC]+", NewScheme())
Dict{String, AbstractScheme} with 1 entry:
  "[+NSC]+" => NewScheme

julia> parse_adduct("[2M+NSC]+")
(adduct = NewScheme, ncore = 2)
```
"""
set_scheme!(nm::AbstractString, scheme::AbstractScheme) = push!(scheme_name(), nm => scheme)