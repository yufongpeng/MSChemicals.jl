# Interface AbstractChemical
@info "Defining interfaces for custom chemical and scheme types"

abstract type Hexose <: AbstractChemical end
struct Glucose <: Hexose
    chirality::String
    nD::Int
    n13C::Int
    retentiontime::Float64
end
struct Galactose <: Hexose
    chirality::String
    nD::Int
    n13C::Int
    retentiontime::Float64
end
struct Mannose <: Hexose 
    chirality::String
    nD::Int
    n13C::Int
    retentiontime::Float64
end
copy(x::Glucose) = Glucose(x.chirality, x.nD, x.n13C, x.retentiontime)
copy(x::Galactose) = Galactose(x.chirality, x.nD, x.n13C, x.retentiontime)
copy(x::Mannose) = Mannose(x.chirality, x.nD, x.n13C, x.retentiontime)
function hash(x::Hexose, h::UInt) 
    h = hash(x.chirality, h) 
    h = hash(x.nD, h) 
    h = hash(x.n13C, h) 
    hash(x.retentiontime, h) 
end
==(x::T, y::T) where {T <: Hexose} = x.chirality == y.chirality && x.nD == y.nD && x.n13C == y.n13C && x.retentiontime == y.retentiontime


repr_isotope(m::Hexose) = repr_isotope(m.nD, m.n13C)
function repr_isotope(nD, n13C)
    if nD == 0 && n13C == 0
        ""
    elseif nD == 0
        string("[13C", n13C, "]")
    elseif n13C == 0
        string("[D", nD, "]")
    else
        string("[13C", n13C, ",D", nD, "]")
    end
end

chemicalname(m::T; kwargs...) where {T <: Hexose} = string(m.chirality, "-", string(T), repr_isotope(m))
chemicalelements(m::Hexose; kwargs...) = 
    [
        "C"     => 6 - m.n13C,
        "[13C]" => m.n13C,
        "H"     => 12 - m.nD, 
        "D"     => m.nD,
        "O"     => 6
    ]

chemicalabbr(m::Hexose; kwargs...) = string(first(chemicalname(m), 5), repr_isotope(m))
chemicalabbr(m::Glucose; kwargs...) = string(first(chemicalname(m), 4), "c", repr_isotope(m))
chemicalsmiles(m::Hexose; kwargs...) = ""

# Fatty acid 
struct FattyAcid <: AbstractChemical
    ncb::Int 
    ndb::Int
    nD::Int 
    n13C::Int 
end

# PS     
struct DiacylPS <: AbstractChemical
    headgroup::Tuple{Int, Int}
    fa1::FattyAcid
    fa2::FattyAcid
    retentiontime::Float64
end
copy(x::FattyAcid) = FattyAcid(x.ncb, x.ndb, x.nD, x.n13C)
copy(x::DiacylPS) = DiacylPS(x.headgroup, copy(x.fa1), copy(x.fa2), x.retentiontime)
function hash(x::FattyAcid, h::UInt) 
    h = hash(x.ncb, h) 
    h = hash(x.ndb, h) 
    h = hash(x.nD, h) 
    hash(x.n13C, h) 
end
function hash(x::DiacylPS, h::UInt) 
    h = hash(x.headgroup, h) 
    h = hash(x.fa1, h) 
    h = hash(x.fa2, h) 
    hash(x.retentiontime, h) 
end
==(x::FattyAcid, y::FattyAcid) = x.ncb == y.ncb && x.ndb == y.ndb && x.nD == y.nD && x.n13C == y.n13C 
==(x::DiacylPS, y::DiacylPS) = x.headgroup == y.headgroup && x.fa1 == y.fa1 && x.fa2 == y.fa2 && x.retentiontime == y.retentiontime

repr_headgroup(m::DiacylPS) = string("PS", repr_isotope(m.headgroup...))
repr_fa(fa::FattyAcid) = string(fa.ncb, ":", fa.ndb, repr_isotope(fa.nD, fa.n13C))
chemicalname(m::FattyAcid; kwargs...) = string("FA", " ", repr_fa(m))
chemicalname(m::DiacylPS; kwargs...) = string(repr_headgroup(m), " ", repr_fa(m.fa1), "/", repr_fa(m.fa2))
function chemicalelements(m::FattyAcid; kwargs...)
    filter(x -> last(x) > 0, [
        "C"     => chemicalncb(m) - chemicaln13C(m), 
        "[13C]" => chemicaln13C(m),
        "H"     => chemicalncb(m) * 2 - chemicalndb(m) * 2 - chemicalnD(m),
        "D"     => chemicalnD(m),
        "O"     => 2
    ])
end
chemicaln13C(m::FattyAcid; kwargs...) = getchemicalproperty(m, :n13C, 0)
chemicalncb(m::FattyAcid; kwargs...) = getchemicalproperty(m, :ncb, 0)
chemicalnD(m::FattyAcid; kwargs...) = getchemicalproperty(m, :nD, 0)
chemicalndb(m::FattyAcid; kwargs...) = getchemicalproperty(m, :ndb, 0)

function chemicalelements(m::DiacylPS; kwargs...)
    fa1 = chemicalelements(m.fa1)
    fa2 = chemicalelements(m.fa2)
    i = findfirst(x -> ==(first(x), "H"), fa1)
    fa1[i] = "H" => (last(fa1[i]) - 1)
    i = findfirst(x -> ==(first(x), "H"), fa2)
    fa2[i] = "H" => (last(fa2[i]) - 1)
    i = findfirst(x -> ==(first(x), "O"), fa1)
    fa1[i] = "O" => (last(fa1[i]) - 1)
    i = findfirst(x -> ==(first(x), "O"), fa2)
    fa2[i] = "O" => (last(fa2[i]) - 1)
    filter(x -> last(x) > 0, vcat(["C" => (6 - m.headgroup[2]), "[13C]" => m.headgroup[2], "H" => (12 - m.headgroup[1]), "D" => m.headgroup[1], "N" => 1, "O" => 8, "P" => 1], fa1, fa2))
end
chemicalformula(m::DiacylPS; unique = false, kwargs...) = chemicalformula(chemicalelements(m); unique, kwargs...)
chemicalsmiles(m::FattyAcid; kwargs...) = ""
chemicalsmiles(m::DiacylPS; kwargs...) = ""
chemicalnD(m::DiacylPS; kwargs...) = +(m.headgroup[1], chemicalnD(m.fa1; kwargs...), chemicalnD(m.fa2; kwargs...))
chemicaln13C(m::DiacylPS; kwargs...) = +(m.headgroup[2], chemicaln13C(m.fa1; kwargs...), chemicaln13C(m.fa2; kwargs...))
chemicalncb(m::DiacylPS; kwargs...) = chemicalncb(m.fa1; kwargs...) + chemicalncb(m.fa2; kwargs...)
chemicalndb(m::DiacylPS; kwargs...) = chemicalndb(m.fa1; kwargs...) + chemicalndb(m.fa2; kwargs...)

# Interface AbstractScheme
struct LossSerine <: AbstractStructuralScheme end
struct LossProtonSerine <: AbstractStructuralScheme end
struct Serine{T} <: AbstractChemicalWrapper{T} 
    chemical::T 
end
copy(x::LossSerine) = x 
copy(x::LossProtonSerine) = x 

function Serine(; nD = 0, n13C = 0) 
    nD == 0 && n13C == 0 && return Serine(Chemical("Serine", "C3H5NO2"; abbreviation = "Ser"))
    @assert 0 <= n13C < 4
    nC = 3 - n13C 
    reprC = nC > 1 ? string("C", nC) : nC > 0 ? "C" : ""
    repr13C = n13C > 1 ? string("[13C]", n13C) : n13C > 0 ? "[13C]" : ""
    @assert 0 <= nD < 5
    nH = 5 - nD 
    reprH = nH > 1 ? string("H", nH) : nH > 0 ? "H" : ""
    reprD = nD > 1 ? string("D", nD) : nH > 0 ? "D" : ""
    nD == 0 && return Serine(Chemical("Serine[13C$(n13C > 1 ? n13C : "")]", "$(reprC)$(repr13C)H5NO2"; abbreviation = "Ser[13C$(n13C > 1 ? n13C : "")]"))
    n13C == 0 && return Serine(Chemical("Serine[D$(nD > 1 ? nD : ""))]", "C3$(reprH)$(reprD)NO2"; abbreviation = "Ser[D$(nD > 1 ? nD : ""))]"))
    Serine(Chemical("Serine[D$(nD > 1 ? nD : ""),13C$(n13C > 1 ? n13C : "")]", "$(reprC)$(repr13C)$(reprH)$(reprD)NO2"; abbreviation = "Ser[D$(nD > 1 ? nD : ""),13C$(n13C > 1 ? n13C : "")]"))
end
set_scheme!("[-H-Ser]-", LossProtonSerine())
set_scheme!("[-Ser-H]-", LossProtonSerine())
set_schabbr!("Ser", Serine())

const SES = StructuralElementalScheme
const LossProton = ElementalScheme{false, <:Proton}
const DiacylPSmH = AdductIon{<:DiacylPS, <:SES{<:LossProton}}
const DiacylPSmS = AdductIon{<:DiacylPS, <:SES{LossSerine}}
const DiacylPSmHS = AdductIon{<:DiacylPS, <:SES{LossProtonSerine}}

# completescheme(precursor::DiacylPS, product::LossSerine) = SES(product, ChemicalLoss(Serine(; nD = first(precursor.headgroup), n13C = last(precursor.headgroup))))
# completescheme(precursor::DiacylPSmH, product::LossSerine) = SES(product, ChemicalLoss(Serine(; nD = first(ioncore(precursor).headgroup), n13C = last(ioncore(precursor).headgroup))))
# completescheme(precursor::DiacylPS, product::LossProtonSerine) = SES(product, ChemicalLoss(AdductIon(Serine(; nD = first(precursor.headgroup), n13C = last(precursor.headgroup)), ChemicalGain(Proton()))))

elementalscheme(precursor::DiacylPS, product::LossSerine) = ChemicalLoss(Serine(; nD = first(precursor.headgroup), n13C = last(precursor.headgroup)))
elementalscheme(precursor::DiacylPSmH, product::LossSerine) = ChemicalLoss(Serine(; nD = first(ioncore(precursor).headgroup), n13C = last(ioncore(precursor).headgroup)))
elementalscheme(precursor::DiacylPS, product::LossProtonSerine) = ChemicalLoss(AdductIon(Serine(; nD = first(precursor.headgroup), n13C = last(precursor.headgroup)), ChemicalGain(Proton())))

# `LossSerine()` -> `LossProtonSerine()` for 
# adductionscheme(precursor::DiacylPSmH, product::SES{LossSerine}) = completescheme(ioncore(precursor), LossProtonSerine())
adductionscheme(precursor::DiacylPSmH, product::LossSerine) = completescheme(ioncore(precursor), LossProtonSerine())
# adductionscheme(precursor::DiacylPSmS, product::SES{<:LossProton}) = completescheme(ioncore(precursor), LossProtonSerine())
adductionscheme(precursor::DiacylPSmS, product::LossProton) = completescheme(ioncore(precursor), LossProtonSerine())

chemicalname(::LossSerine; n = 1, loss = false, delim = "", kwargs...) = string(delim, loss ? "Gain_" : "Loss_", n > 1 ? n : "", "Serine")
chemicalname(::LossProtonSerine; n = 1, loss = false, delim = "", kwargs...) = string(delim, loss ? "Gain_" : "Loss_", n > 1 ? n : "", "[Serine+H]+")
chemicalabbr(::LossSerine; n = 1, loss = false, delim = "", kwargs...) = string(delim, loss ? "+" : "-", n > 1 ? n : "", "Ser")
chemicalabbr(::LossProtonSerine; n = 1, loss = false, delim = "", kwargs...) = string(delim, loss ? "+" : "-", n > 1 ? n : "", "[Ser+H]+")

struct SN1Acyl <: StructuralChemicalScheme end
struct SN2Acyl <: StructuralChemicalScheme end
# completescheme(precursor::DiacylPSmH, product::SN1Acyl) = SES(product, AdductIon(ioncore(precursor).fa1, ChemicalLoss(Proton())))
# completescheme(precursor::DiacylPSmH, product::SN2Acyl) = SES(product, AdductIon(ioncore(precursor).fa2, ChemicalLoss(Proton())))
# completescheme(precursor::DiacylPSmHS, product::SN1Acyl) = SES(product, AdductIon(ioncore(precursor).fa1, ChemicalLoss(Proton())))
# completescheme(precursor::DiacylPSmHS, product::SN2Acyl) = SES(product, AdductIon(ioncore(precursor).fa2, ChemicalLoss(Proton())))
# detectedchemical(precursor::DiacylPSmH, product::SN1Acyl) = AdductIon(ioncore(precursor).fa1, ChemicalLoss(Proton()))
# detectedchemical(precursor::DiacylPSmH, product::SN2Acyl) = AdductIon(ioncore(precursor).fa2, ChemicalLoss(Proton()))
# detectedchemical(precursor::DiacylPSmHS, product::SN1Acyl) = AdductIon(ioncore(precursor).fa1, ChemicalLoss(Proton()))
# detectedchemical(precursor::DiacylPSmHS, product::SN2Acyl) = AdductIon(ioncore(precursor).fa2, ChemicalLoss(Proton()))
elementalscheme(precursor::DiacylPSmH, product::SN1Acyl) = AdductIon(ioncore(precursor).fa1, ChemicalLoss(Proton()))
elementalscheme(precursor::DiacylPSmH, product::SN2Acyl) = AdductIon(ioncore(precursor).fa2, ChemicalLoss(Proton()))
elementalscheme(precursor::DiacylPSmHS, product::SN1Acyl) = AdductIon(ioncore(precursor).fa1, ChemicalLoss(Proton()))
elementalscheme(precursor::DiacylPSmHS, product::SN2Acyl) = AdductIon(ioncore(precursor).fa2, ChemicalLoss(Proton()))

chemicalname(::SN1Acyl; n = 1, kwargs...) = string(n > 1 ? n : "", "Sn1_Acyl")
chemicalname(::SN2Acyl; n = 1, kwargs...) = string(n > 1 ? n : "", "Sn2_Acyl")
chemicalabbr(::SN1Acyl; n = 1, kwargs...) = string(n > 1 ? n : "", "Sn1")
chemicalabbr(::SN2Acyl; n = 1, kwargs...) = string(n > 1 ? n : "", "Sn2")

struct FC <: AbstractChemical 
    name::String
    formula::String
end

struct FS <: AbstractElementalScheme 
    name::String
    formula::String
end

struct EC <: AbstractChemical 
    name::String
    elements::Vector{Pair{String, Int}}
end

struct ES <: AbstractElementalScheme 
    name::String
    elements::Vector{Pair{String, Int}}
end

@info "Defining custom chemicals and schema"

glc = Glucose("D", 0, 0, 1.5)
gld = Glucose("D", 6, 0, 1.5)
ps = DiacylPS((0, 0), FattyAcid(18, 0, 0, 0), FattyAcid(20, 4, 0, 0), 7.8)
psi1 = DiacylPS((3, 3), FattyAcid(18, 0, 0, 0), FattyAcid(20, 4, 0, 0), 7.8)
psi2 = DiacylPS((0, 0), FattyAcid(18, 0, 5, 0), FattyAcid(20, 4, 0, 0), 7.8)

@info "Defining generic adduct ions"

iglc = [AdductIon(glc, "[M+H]+")]
igld = [AdductIon(gld, "[M+H]+")]
ips = [AdductIon(ps, "[M-H-Ser]-"), AdductIon(ps, "[2M+H]+"), AdductIon(ps, "[M-H]-")]
ipsi1 = [AdductIon(psi1, "[M-H-Ser]-"), AdductIon(psi1, "[2M+H]+"), AdductIon(psi1, "[M-H]-")]
ipsi2 = [AdductIon(psi2, "[M-H-Ser]-"), AdductIon(psi1, "[2M+H]+"), AdductIon(psi2, "[M-H]-")]

@info "Defining chemical transitions"

sp1 = ChemicalSeries(ipsi1[1], SN2Acyl())
sp2 = ChemicalSeries(ipsi1[3] => LossSerine() => SN2Acyl())
sp3 = ChemicalSeries(ipsi2[1], SN1Acyl())
sp4 = ChemicalSeries(ipsi2[3] => LossSerine() => SN1Acyl())
sp5 = ChemicalSeries(isotopomerize(ipsi1[1], ["[13C]" => 5]) => isotopomerize(completescheme(ipsi1[1], outputchemical(sp1)), ["[13C]" => 5]))