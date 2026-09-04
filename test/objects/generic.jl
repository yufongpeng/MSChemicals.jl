@info "Defining generic chemicals and schema"

cglc = Chemical("Glucose", ["C" => 6, "H" => 12, "O" => 6]; retentiontime = 1.5, abbreviation = "Glc", SMILES = "")
fglc = FormulaChemical(["C" => 6, "H" => 12, "O" => 6]; retentiontime = 1.5, abbreviation = "Glc", SMILES = "")
cgld = Chemical("Glucose_d6", "C6H6D6O6"; retentiontime = 1.5, abbreviation = "Glc[D6]", SMILES = "")
cps = Chemical("PS 18:0/20:4", "C44H78NO10P"; retentiontime = 7.8)
cpsi1 = Chemical("PS[D3,13C3] 18:0/20:4", "C41[13C]3H75D3NO10P"; retentiontime = 7.8)
cpsi2 = Chemical("PS 18:0[D5]/20:4(5Z,8Z,11Z,14Z)", "C44H73D5NO10P"; retentiontime = 7.8)
fa1 =  Chemical("FA 18:0", "C18H36O2")
fa1i =  Chemical("FA 18:0[D5]", "C18D5H31O2")
cserine = Chemical("Serine", "C3H5NO2"; abbreviation = "Ser")
cserinei = Chemical("Serine[D3,13C3]", "[13C]3H2D3NO2"; abbreviation = "Ser[D3,13C3]")
lossserine = ChemicalLoss(cserine)
lossserinei = ChemicalLoss(cserinei)
losshserine = ChemicalLoss(AdductIon(cserine, "[M+H]+"))
losshserinei = ChemicalLoss(AdductIon(cserinei, "[M+H]+"))
losshserinegainwater = ChemicalSchema(ChemicalLoss(AdductIon(cserine, "[M+H]+")), ChemicalGain(Water()))
losshserinegaincolossco = ChemicalSchema(ChemicalLoss(AdductIon(cserine, "[M+H]+")), parse_chemical("+CO"), parse_chemical("-CO"))
# Generic structure interface
# Mix `[M-H]-` and neutral loss `lossserine` into `losshserine` for `AddductIon(cps, "[M-H]-")`
push!(cps.property, :schema => [
    ChemicalLoss(Proton()) => [
        lossserine => losshserine
    ] 
])
# Change neutral loss `lossserine` to `lossserinei` for `cpsi1`
# Change neutral loss `losshserine` to `losshserinei` for `cpsi1`
# Change neutral loss `lossserine` to `lossserinei` for `AddductIon(cpsi1, "[M-H]-")`
push!(cpsi1.property, :structure => [
    nothing => [
        lossserine => lossserinei,
        losshserine => losshserinei
    ],
    ChemicalLoss(Proton()) => [
        lossserine => lossserinei  
    ]
])
# Mix `[M-H]-` and neutral loss `lossserine` into `losshserine` for `AddductIon(cpsi1, "[M-H]-")`
push!(cpsi1.property, :schema => [
    ChemicalLoss(Proton()) => [
        lossserine => losshserine
    ] 
])
# Change product `AdductIon(fa1, "[M-H]-")` to `AdductIon(fa1i, "[M-H]-")` for `AddductIon(cpsi2, "[M-H]-")`
push!(fa1.property, :chemicalscheme => [
    ChemicalLoss(Proton()) => :sn1fa
])
push!(fa1i.property, :chemicalscheme => [
    ChemicalLoss(Proton()) => :sn1fa
])
push!(cpsi2.property, :structure => [
    ChemicalLoss(Proton()) => [
        :sn1fa => AdductIon(fa1i, ChemicalLoss(Proton()))     
    ],
    losshserine => [
        :sn1fa => AdductIon(fa1i, ChemicalLoss(Proton()))     
    ]
])
# Mix `[M-H]-` and neutral loss `lossserine` into `losshserine` for `AddductIon(cpsi2, "[M-H]-")`
push!(cpsi2.property, :schema => [
    ChemicalLoss(Proton()) => [
        lossserine => losshserine
    ] 
])

@info "Defining generic adduct ions"

icglc = [AdductIon(cglc, "[M+H]+"), AdductIon(cglc, "[M+H-H2O]+")]
icgld = [AdductIon(cgld, "[M+H]+"), AdductIon(cgld, "[M+H-H2O]+")]
icps = [ionize(cps; adduct = losshserine), AdductIon(cps, "[2M+H]+"), AdductIon(cps, "[M-H]-")]
icpsi1 = [ionize(cpsi1, losshserine), AdductIon(cpsi1, "[2M+H]+"), AdductIon(cpsi1, "[M-H]-")]
icpsi2 = [ionize(cpsi2; adduct = losshserine), AdductIon(cpsi2, "[2M+H]+"), AdductIon(cpsi2, "[M-H]-")]

@info "Defining chemical transitions"

cp0 = ChemicalSeries(icps[3] => lossserine)
cp1 = ChemicalSeries(icps[1], AdductIon(fa1, ChemicalLoss(Proton())))
cp2 = ChemicalSeries(icps[3] => lossserine => outputchemical(cp1))
cp3 = ChemicalSeries(icps[3], cp1)
cp4 = ChemicalSeries(ionize(isotopomerize(icps[3], ["[13C]" => 5]); adduct = ElementalScheme(false, Proton())) => isotopomerize(completescheme(icps[1], outputchemical(cp1)), ["[13C]" => 5]))
cp5 = ChemicalSeries(ionize(isotopomerize(isotopomerize(icpsi1[3], ["[13C]" => 5]), ["D" => 2]), ElementalScheme(false, Proton())), ChemicalSeries(isotopomerize(completescheme(icpsi1[3], lossserine), ["D" => 1]), isotopomerize(completescheme(icpsi1[1], outputchemical(cp1)), ["[13C]" => 5])))
cp6 = ChemicalSeries(icpsi2[3] => ionize(fa1, "[M-H]-"))
cp7 = ChemicalSeries(icpsi1[3], lossserine) 
pt1 = peak_table(MSScan(Isotopologues(icglc[1]; abundance = 1e5, threshold = crit(1e1, 1e-2))))

@info "Defining criteria"

ct1 = crit(10)
ct2 = rcrit(0.2)
ct3 = crit(10, 0.2)
qualified_peak1(x, x̂, ct) = all(c -> in(x, c), MSC.makecrit_delta(ct, x̂))
qualified_peak2(x, x̂, ct) = any(c -> x >= c, MSC.makecrit_value(ct, x̂))
