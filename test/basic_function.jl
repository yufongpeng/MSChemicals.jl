@info "Running tests for Base interface"

@testset "Base interface" begin
    @test copy(cglc) == cglc
    @test hash(copy(cglc)) == hash(cglc) 
    @test copy(fglc) == fglc
    @test hash(copy(fglc)) == hash(fglc) 
    @test copy(lossserine) == ChemicalLoss(cserine)
    @test hash(copy(lossserine)) == hash(lossserine)
    @test copy(structuralscheme(ionadduct(icpsi1[2]))) == structuralscheme(ionadduct(icpsi1[2]))
    @test hash(copy(structuralscheme(ionadduct(icpsi1[2])))) == hash(structuralscheme(ionadduct(icpsi1[2])))
    @test copy(last(icglc)) == last(icglc)
    @test hash(copy(last(icglc))) == hash(last(icglc))
    @test copy(cp3) == cp3
    @test hash(copy(cp3)) == hash(cp3)
    @test copy(pt1.Chemical[3]) == pt1.Chemical[3]
    @test hash(copy(pt1.Chemical[3])) == hash(pt1.Chemical[3])
    @test copy(it3.Chemical[2]) == it3.Chemical[2]
    @test hash(copy(git3.Chemical[2])) == hash(git3.Chemical[2])
    @test copy(git3.Chemical[2]) == git3.Chemical[2]
    @test hash(copy(git3.Chemical[2])) == hash(git3.Chemical[2])
    @test copy(chemicaltransition(itit12.Chemical[11])[2]) == chemicaltransition(itit12.Chemical[11])[2]
    @test hash(copy(chemicaltransition(itit12.Chemical[11])[2])) == hash(chemicaltransition(itit12.Chemical[11])[2])
    @test copy(chemicaltransition(gitit12.Chemical[7])[2]) == chemicaltransition(gitit12.Chemical[7])[2]
    @test hash(copy(chemicaltransition(gitit12.Chemical[7])[2])) == hash(chemicaltransition(gitit12.Chemical[7])[2])
end

@info "Running tests for Default chemical parsing"

@testset "Default chemical parsing" begin
    @test ischemicalequal(parse_chemical(cglc), cglc)
    @test !(@test_noerror parse_chemical(ChemicalSchemeParser(), "C6H12O"))
    @test ischemicalequal(parse_chemical(ChemicalEntityParser(ChemicalParser()), "[Glucose_d6+H]+"; formula = "C6H6D6O6", retentiontime = 1.5), AdductIon(parse_chemical(ChemicalParser(; formula = "C6H12O6"), "Glucose_d6"; retentiontime = 1.5, formula = "C6H6D6O6"), "[M+H]+"))
    @test ischemicalequal(parse_chemical("[+Na+H-H2O]2+"), parse_chemical("[-H2O+Na+H]2+"))
    @test ischemicalequal(parse_chemical("-[2H]2+"), parse_chemical("[-2H]2-"))
    @test ischemicalequal(parse_chemical("+[NH4]+"), parse_chemical("[+NH4]+"))
    @test ischemicalequal(parse_chemical("C6H12O6"), FormulaChemical("C6H12O6"))
    @test ischemicalequal(parse_chemical(ChemicalTransitionParser(ChemicalEntityParser(FormulaChemicalParser(; retentiontime = 5))), "[C6H12O6]+"; chemicaltype = "Organic", retentiontime = 1.5), AdductIon(FormulaChemical("C6H12O6"; chemicaltype = "Organic", retentiontime = 1.5), "[M]+"))
    @test ischemicalequal(parse_chemical("[C6H12O6+H]+"), AdductIon(FormulaChemical("C6H12O6"), ChemicalGain(Proton())))
    @test ischemicalequal(parse_chemical("[C6H12O6+H]+" => "-H2O"), parse_chemical(AdductIon(FormulaChemical("C6H12O6"), ChemicalGain(Proton())), ChemicalLoss(Water())))
    @test ischemicalequal(detectedchemical(parse_chemical("[C6H12O6+H]+" => "-H2O")), parse_chemical("[C6H12O6+H-H2O]+"))
    @test ischemicalequal(parse_chemical(ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1)), "[C6H12O6+2H]2+", "C6H10O5"), parse_chemical(AdductIon(FormulaChemical("C6H12O6"), "[M+2H]2+") => AdductIon(FormulaChemical("C6H10O5"), "[M]+")))
    @test ischemicalequal(parse_chemical(ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 1)), "[C6H12O6+2H]2+ -> -H3O ->  -[CO2]"), parse_chemical(ChemicalSeries(AdductIon(FormulaChemical("C6H12O6"), "[M+2H]2+"), ChemicalSeries(ChemicalLoss(AdductIon(FormulaChemical("H3O"), "[M]+")), ChemicalLoss(FormulaChemical("CO2"))))))
    @test ischemicalequal(parse_chemical(ChemicalTransitionParser(ChemicalExpressionParser(; charge = 1, loss = 1)), "[C6H12O6+2H]2+" => "-H3O" => "-[CO2]"), parse_chemical(AdductIon(FormulaChemical("C6H12O6"), "[M+2H]2+") => ChemicalSeries(ChemicalLoss(AdductIon(FormulaChemical("H3O"), "[M]+")), ChemicalLoss(FormulaChemical("CO2")))))
    @test parse_adduct((1, 2)) == (1, 2)
    @test parse_adduct(AdductParser(true), (1, 2)) == parse_adduct(AdductParser(true), (adduct = 1, ncore = 2))
    @test parse_adduct((adduct = 1, ncore = 2)) == parse_adduct(AdductParser(false), (adduct = 1, ncore = 2))
end