@testset "Attributes" begin
    @testset "Basic chemicalformula and elements" begin 
        @test MSC.gain_elements(MSC.dictionary_elements(MSC.ElementsVector(["D"], [5])), Dictionary(["H"], [5])) == MSC.loss_elements(Dict("D" => 5, "H" => 10), MSC.dictionary_elements(Dictionary, MSC.ElementsVector(["H"], [5])))
        @test MSC.gain_elements(["D" => 5], MSC.dictionary_elements(Dictionary, Dictionary(["H"], [5]))) == MSC.unique_elements(MSC.loss_elements(["D" => 5, "H" => 10], MSC.unique_elements(Dictionary, Dict("H" => 5))))
        @test MSC.gain_elements(Dictionary(["D"], [5]), MSC.dictionary_elements(Dict, Dictionary(["H"], [5]))) == MSC.loss_elements(MSC.unique_elements(Dictionary, ["D" => 5, "H" => 10]), MSC.dictionary_elements(Dict, Dict("H" => 5)))
        @test MSC.gain_elements(MSC.unique_elements(Vector{Pair}, Dictionary(["D"], [5])), MSC.unique_elements(Dict, Dictionary(["H"], [5]))) == MSC.unique_elements(Vector{Pair}, MSC.loss_elements(MSC.unique_elements(Dictionary, Dictionary(["D", "H"], [5, 10])), MSC.unique_elements(Vector{Pair}, Dict("H" => 5))))
        @test chemicalformula(["C" => 2, "H" => 5, "O" => 1, "H" => 1]; unique = false) == "C2H5OH"
        @test chemicalformula(["C" => 2, "H" => 5, "O" => 1, "H" => 1]; unique = false, ischemical = false) == "+C2H5OH"
        @test chemicalformula(["C" => -2, "H" => -5, "O" => -1, "H" => -1]; unique = false, ischemical = false, loss = false) == "-C2H5OH"
        @test chemicalformula(["C" => 2, "H" => 5, "O" => 1, "H" => 1, "C" => -1, "O" => -1]; unique = false, loss = false) == "CH6"
        @test chemicalformula(["C" => -2, "H" => -5, "O" => -1, "H" => -1, "C" => 1, "O" => 1, "F" => 1]; unique = false, ischemical = false, loss = false) == "-C2H5OH+COF"
        @test chemicalformula(MSC.dictionary_elements(Dictionary, ["C" => -2, "H" => -5, "O" => -1, "H" => -1, "C" => 1, "O" => 1, "F" => 1]); ischemical = false, loss = false) == "-CH6+F"
        @test chemicalformula(EC("C6H12", ["C" =>6, "H" => 12])) == chemicalformula(FC("C6H12", "C6H12"))
        @test chemicalelements(EC("C6H12", ["C" => 6, "H" => 12])) == chemicalelements(FC("C6H12", "C6H12"))
        @test chemicalformula(ES("-C6H12", ["C" => -6, "H" => -12]); loss = true) == chemicalformula(FS("C6H12", "C6H12"))
        @test chemicalformula(ES("-C6H12", ["C" => -6, "H" => -12])) == chemicalformula(FS("C6H12", "C6H12"); loss = true)
        @test chemicalelements(ES("-C6H12", ["C" => -6, "H" => -12])) == chemicalelements(FS("C6H12", "C6H12"); loss = true)
        @test chemicalname(ES("-C6H12", ["C" => -6, "H" => -12])) == chemicalname(FS("-C6H12", "-C6H12"))
    end
    @testset "getchemicalproperty" begin
        @test getchemicalproperty(cglc, :name, "MISSING") == "Glucose"
        @test getchemicalproperty(cglc, :doesnotexist, "DEFAULT") == "DEFAULT"
        @test getchemicalproperty(icglc[1], :charge, 0) == charge(ioncore(icglc[1]))
        @test getchemicalproperty(icglc[1], :abbreviation, "DEFAULT") == chemicalabbr(ioncore(icglc[1]))
    end

    @testset "Basic chemical attributes" begin
        @test @test_noerror test_show(cglc)

        @test chemicalname(cglc) == "Glucose"
        @test chemicalabbr(cglc) == "Glc"
        @test chemicalformula(cglc) == "C6H12O6"
        @test MSC.dictionary_elements(chemicalelements(cglc)) == Dict("C" => 6, "H" => 12, "O" => 6)

        @test chemicalname(cgld) == "Glucose_d6"
        @test chemicalabbr(cgld) == "Glc[D6]"
        @test chemicalformula(cgld) == "C6H6D6O6"
        @test MSC.dictionary_elements(chemicalelements(cgld)) == Dict("C" => 6, "D" => 6, "H" => 6, "O" => 6)

        @test chemicalname(fglc) == "C6H12O6"
        @test chemicalabbr(fglc) == "Glc"
        @test chemicalformula(fglc) == "C6H12O6"
        @test MSC.dictionary_elements(chemicalelements(fglc)) == Dict("C" => 6, "H" => 12, "O" => 6)

        @test chemicalname(cps) == "PS 18:0/20:4"
        @test chemicalabbr(cps) == chemicalname(cps)
        @test chemicalformula(cps) == "C44H78NO10P"
        @test MSC.dictionary_elements(chemicalelements(cps)) == Dict("C" => 44, "H" => 78, "N" => 1, "O" => 10, "P" => 1)

        @test isnan(mz(cglc))
        @test isapprox(mmi(cpsi1) * 2 + mmi("H") - MSC.ME, mz(cpsi1, "[2M+H]+"))
        @test isapprox(molarmass(cpsi1), molarmass(chemicalformula(cpsi1)))
    end

    @testset "Adduct ion and wrapper attributes" begin
        @test ioncore(icglc[1]) == cglc
        @test !isnothing(ionadduct(icglc[1]))
        @test isapprox(mz(icpsi1[1]), mz(icps[1]))
        @test ncore(icpsi1[2]) == 2
        @test charge(icglc[1]) == 1
        @test ncharge(cglc) == 0
        @test isempty(chemicalsmiles(icglc[1]))
        @test retentiontime(icglc[1]) == retentiontime(cglc)
        @test chemicalformula(icglc[1]) == chemicalformula(chemicalelements(icglc[1]); ischemical = true)
        @test chemicalformula(icglc[1]; loss = false) == chemicalformula(chemicalelements(icglc[1]; loss = false); ischemical = true)
        @test chemicalname(icglc[1]; bracket = false) == string(chemicalname(cglc; bracket = false), chemicalabbr(ionadduct(icglc[1]); bracket = false))
        @test ischemicalequal(icpsi1[1], AdductIon(cpsi1, losshserine))
        @test isapprox(retentiontime(icpsi2[1]), 7.8) 
        # chemicalstructure  
        @test ischemicalequal(chemicalentity(icgld[1]), icgld[1])
        @test ischemicalequal(first(chemicalspecies(icgld[1])), icgld[1])
        @test ischemicalequal(first(chemicaltransition(icgld[1])), icgld[1])
        @test ischemicalequal(chemicalparent(icgld[1]), icgld[1])
        @test isempty(isotopomersisotopes(icgld[1]))
        @test ischemicalequal(analyzedchemical(icgld[1]), icgld[1])
        @test ischemicalequal(detectedchemical(icgld[1]), icgld[1])
        @test ischemicalequal(inputchemical(icgld[1]), icgld[1])
        @test ischemicalequal(outputchemical(icgld[1]), icgld[1])
        @test ischemicalequal(first(seriesanalyzedchemical(icgld[1])), first(chemicaltransition(icgld[1])))
        @test isempty(detectedisotopes(icgld[1]))
        @test isempty(last(seriesanalyzedisotopes(icgld[1])))
        @test detectedcharge(icgld[1]) == last(seriesanalyzedcharge(icgld[1]))
        @test detectedelements(icgld[1]) == last(seriesanalyzedelements(icgld[1]))
    end

    @testset "Scheme and loss/gain attributes" begin
        @test chemicalname(lossserine) == "Loss_Serine"
        @test chemicalabbr(lossserine) == "-Ser"
        @test chemicalformula(lossserine) == chemicalformula(chemicalelements(lossserine); ischemical = false)
        @test MSC.dictionary_elements(chemicalelements(lossserine)) == Dict("C" => -3, "H" => -5, "N" => -1, "O" => -2)
        @test isempty(isotopomersisotopes(lossserine))

        @test chemicalname(losshserine) == "Loss_[Serine+H]+"
        @test chemicalabbr(losshserine) == "-[Ser+H]+"
        @test chemicalformula(losshserine) == chemicalformula(chemicalelements(losshserine); ischemical = false)
        @test isempty(isotopomersisotopes(losshserine))

        @test chemicalname(losshserinegainwater) == "Loss_[Serine+H]+|Gain_Water"
        @test chemicalabbr(losshserinegainwater) == "[-Ser-H+H2O]-"
        @test chemicalformula(losshserinegainwater) == chemicalformula(chemicalelements(losshserinegainwater); ischemical = false)
        @test isempty(isotopomersisotopes(losshserinegainwater))

        @test chemicalname(lossserine; loss = true) == "Gain_Serine"
        @test chemicalabbr(lossserine; loss = true) == "+Ser"

        @test chemicalname(losshserinegaincolossco) == "Loss_[Serine+H]+|+[CO]-[CO]"
        @test chemicalabbr(losshserinegaincolossco) == "[-Ser-H+CO-CO]-"
        @test chemicalformula(losshserinegaincolossco) == chemicalformula(chemicalelements(losshserinegaincolossco); ischemical = false)
        @test isempty(isotopomersisotopes(losshserinegaincolossco))

        @test isapprox(mmi(ionadduct(icpsi1[1])), mz(ionadduct(icpsi1[1])))
        @test isapprox(molarmass(ionadduct(icpsi1[1])), molarmass(losshserinei))
        @test chemicalparent(losshserinegaincolossco) == losshserinegaincolossco
        @test elementalscheme(losshserinegaincolossco) == losshserinegaincolossco
        @test structuralscheme(losshserinegaincolossco) == losshserinegaincolossco
    end

    @testset "Custom chemical types" begin
        @test chemicalname(glc) == "D-Glucose"
        @test MSC.unique_elements(chemicalelements(glc)) == ["C" => 6, "H" => 12, "O" => 6]
        @test chemicalformula(glc; unique = true) == "C6H12O6"
        @test chemicalsmiles(glc) == ""

        @test chemicalname(iglc[1]) == "[D-Glucose+H]+"
        @test chemicalabbr(iglc[1]) == "[D-Glc+H]+"
        @test charge(iglc[1]) == 1
        @test ncore(iglc[1]) == 1

        @test chemicalname(ps) == "PS 18:0/20:4"
        @test chemicalformula(ps) == "C6H12NO8PC18H35OC20H31O"
        @test chemicalname(ips[1]) == "[(PS 18:0/20:4)-Ser-H]-"
        @test chemicalname(ionadduct(first(ipsi1))) == "Loss_[Serine[D3,13C3]+H]+"
        @test chemicalabbr(ionadduct(first(ipsi1))) == "-[Ser[D3,13C3]+H]+"
        @test charge(ips[1]) == -1

        @test chemicalentity(Serine()) == Serine()
        @test isnan(retentiontime(Serine()))
        @test isempty(chemicalsmiles(Serine()))
        @test isempty(isotopomersisotopes(Serine()))

        @test parse_chemical("[-Ser-H]-") == LossProtonSerine()
        @test parse_chemical("[-Ser-H+H2O]-") == ChemicalSchema(ChemicalLoss(Serine()), ChemicalLoss(Proton()), ChemicalGain(Water()))
    end

    @testset "Structural scheme attributes" begin
        @test chemicalname(LossProtonSerine()) == "Loss_[Serine+H]+"
        @test chemicalabbr(LossProtonSerine()) == "-[Ser+H]+"
        @test chemicalname(SN1Acyl()) == "Sn1_Acyl"
        @test chemicalabbr(SN1Acyl()) == "Sn1"
        @test chemicalname(SN2Acyl()) == "Sn2_Acyl"
        @test chemicalabbr(SN2Acyl()) == "Sn2"

        @test !(@test_noerror chemicalformula(SN1Acyl(); ischemical = false))
        @test !(@test_noerror chemicalelements(SN1Acyl(); loss = false))
        @test !(@test_noerror charge(SN1Acyl()))
        @test !(@test_noerror elementalscheme(SN1Acyl()))
        @test !(@test_noerror chemicalparent(SN1Acyl()))
        @test !(@test_noerror isotopomersisotopes(SN1Acyl()))
        @test !(@test_noerror groupedisotopomersisotopes(SN1Acyl()))
        @test !(@test_noerror groupedisotopomersabundance(SN1Acyl()))
        @test !(@test_noerror analyzedchemical(SN1Acyl()))
        @test !(@test_noerror mz(SN1Acyl()))
        @test !(@test_noerror molarmass(SN1Acyl()))
    end

    @testset "internal interface of elementalscheme, completescheme, adductionscheme, and detectedchemical" begin 
        @test MSC.elementalscheme(nothing, cglc) == MSC.completescheme(nothing, cglc)
        @test MSC.elementalscheme(cglc, ionadduct(icglc[1])) == MSC.elementalscheme(icglc[1], ionadduct(icglc[1]))
        @test MSC.elementalscheme(cglc, ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water()))) == MSC.elementalscheme(icglc[1], ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water())))
        @test MSC.elementalscheme(cglc, isotopomerize(ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water())), ["[18O]" => 1])) == MSC.elementalscheme(icglc[1], isotopomerize(ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water())), ["[18O]" => 1]))
        @test MSC.elementalscheme(cglc, fa1) == MSC.elementalscheme(icglc[1], fa1)
        @test MSC.elementalscheme(cglc, AdductIon(fa1, "[M-H]-")) == MSC.elementalscheme(icglc[1], AdductIon(fa1, "[M-H]-"))
        @test MSC.elementalscheme(glc, cglc) == MSC.completescheme(glc, cglc)
        @test MSC.elementalscheme(glc, ionadduct(icglc[1])) == MSC.elementalscheme(iglc[1], ionadduct(icglc[1]))
        @test MSC.elementalscheme(glc, ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water()))) == MSC.elementalscheme(iglc[1], ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water())))
        @test MSC.elementalscheme(glc, isotopomerize(ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water())), ["[18O]" => 1])) == MSC.elementalscheme(iglc[1], isotopomerize(ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water())), ["[18O]" => 1]))
        @test MSC.elementalscheme(glc, fa1) == MSC.completescheme(iglc[1], fa1)
        @test MSC.elementalscheme(glc, AdductIon(fa1, "[M-H]-")) == MSC.completescheme(iglc[1], AdductIon(fa1, "[M-H]-"))
        @test MSC.elementalscheme(glc, iglc[1]) == MSC.completescheme(glc, iglc[1])
        @test MSC.completescheme(cglc, cglc) == MSC.completescheme(icglc[1], cglc)
        @test MSC.completescheme(cglc, icglc[1]) == MSC.completescheme(icglc[1], icglc[1])
        @test MSC.completescheme(glc, ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water()))) == chemicalparent(MSC.completescheme(glc, isotopomerize(ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water())), ["[18O]" => 1])))
        @test MSC.completescheme(cglc, ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water()))) == chemicalparent(MSC.completescheme(cglc, isotopomerize(ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water())), ["[18O]" => 1])))
        @test MSC.completescheme(icglc[1], ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water()))) == chemicalparent(MSC.completescheme(icglc[1], isotopomerize(ChemicalSchema(ElementalScheme(true, Proton()), ElementalScheme(false, Water())), ["[18O]" => 1])))
        @test detectedchemical(ps, LossProtonSerine()) == detectedchemical(ps, completescheme(ps, LossProtonSerine()))
        @test detectedchemical(ps, StructuralElementalScheme(SN1Acyl(), ps.fa1)) == detectedchemical(nothing, StructuralElementalScheme(SN1Acyl(), ps.fa1))
        @test !(@test_noerror detectedchemical(ps, SN1Acyl()))
        @test !(@test_noerror detectedchemical(cps, SN1Acyl()))
        @test !(@test_noerror detectedchemical(icps[1], SN1Acyl()))

        @test MSC.adductionscheme(ips[3], ElementalScheme(false, Water())) == ChemicalSchema(ionadduct(ips[3]), StructuralElementalScheme(ElementalScheme(false, Water()), ElementalScheme(false, Water())))
    end

    @testset "Isotopomers" begin
        @test Isotopomers(icglc[1], "C5H12O6[13C]") == it1.Chemical[2]
        @test @test_noerror test_show(it3.Chemical[2])
        @test chemicalname(it3.Chemical[2]) == string(chemicalname(chemicalparent(it3.Chemical[2])), "[13C]")
        @test chemicalabbr(it3.Chemical[2]) == string(chemicalabbr(chemicalparent(it3.Chemical[2])), "[13C]")
        @test chemicalname(itit14.Chemical[2]; bracket = true) == chemicalabbr(itit14.Chemical[2])
        @test chemicalname(itit13.Chemical[4]; bracket = true) == chemicalabbr(itit13.Chemical[4])
        @test chemicalformula(it3.Chemical[2]) == chemicalformula(chemicalelements(it3.Chemical[2]); ischemical = true)
        @test MSC.gain_elements(["D" => 5], isotopomersisotopes(it3.Chemical[2])) == filter(x -> !iselement(first(x)) && last(x) != 0, chemicalelements(it3.Chemical[2]))
        @test chemicalparent(it3.Chemical[2]) == it3.Chemical[2].parent
        @test all(isapprox.(mmi.(it3.Chemical), it3.MZ1; rtol = 20e-6))
        @test all(isapprox.(mz.(it3.Chemical), it3.MZ1; rtol = 20e-6))
        @test charge(it3.Chemical[2]) == -1
        @test isempty(chemicalsmiles(it3.Chemical[2]))
    end

    @testset "Groupedisotopomers" begin
        @test @test_noerror test_show(git3.Chemical[2])
        @test chemicalname(git3.Chemical[2]) == string(chemicalname(chemicalparent(git3.Chemical[2])), "(+1)")
        @test chemicalabbr(git3.Chemical[2]) == chemicalname(git3.Chemical[2])
        @test chemicalname(gitit14.Chemical[2]; bracket = true) == chemicalabbr(gitit14.Chemical[2])
        @test chemicalname(gitit13.Chemical[4]; bracket = true) == chemicalabbr(gitit13.Chemical[4])
        @test chemicalformula(first(chemicalspecies(git3.Chemical[2]))) == chemicalformula(chemicalelements(first(chemicalspecies(git3.Chemical[2]))))
        @test MSC.gain_elements(["D" => 5], isotopomersisotopes(first(chemicalspecies(git3.Chemical[2])))) == filter(!iselement ∘ first, chemicalelements(first(chemicalspecies(git3.Chemical[2]))))
        @test groupedisotopomersabundance(git3.Chemical[2]) == git3.Chemical[2].abundance
        @test chemicalparent(git3.Chemical[2]) == git3.Chemical[2].parent
        @test all(isapprox.(mmi.(git3.Chemical), git3.MZ1; rtol = 20e-6))
        @test all(isapprox.(mz.(git3.Chemical), git3.MZ1; rtol = 20e-6))
        @test ischemicalequal(git3.Chemical[2], git3.Chemical[2])
        @test charge(git3.Chemical[2]) == -1
        @test isempty(chemicalsmiles(git3.Chemical[2]))
    end

    @testset "Round-trip consistency" begin
        for obj in (cglc, cgld, fglc, cps, lossserine, losshserinegainwater, icglc[1], iglc[1], losshserine)
            if obj isa AbstractScheme
                @test chemicalformula(obj) == chemicalformula(chemicalelements(obj); ischemical = false)
            else
                @test chemicalformula(obj) == chemicalformula(chemicalelements(obj); ischemical = true)
            end
        end
    end

    @testset "Isobars" begin
        @test chemicalname(pt1.Chemical[3]) == "Isobars[[Glucose+H]+[18O], [Glucose+H]+[13C2], [Glucose+H]+[17O,13C]]" || chemicalname(pt1.Chemical[3]) == "Isobars[[Glucose+H]+[18O], [Glucose+H]+[13C2], [Glucose+H]+[13C,17O]]"
        @test chemicalabbr(pt1.Chemical[1]) == "Isobars[[Glc+H]+]"
        @test chemicalformula(pt1.Chemical[1]) == chemicalformula(chemicalelements(pt1.Chemical[1]))
        @test isapprox(mz(pt1.Chemical[2], "[M+H]+"), mz(pt1.Chemical[2]))
        @test isapprox(molarmass(pt1.Chemical[1]), molarmass(detectedchemical(pt1.Chemical[1])))
        @test !ischemicalequal(pt1.Chemical[1], pt1.Chemical[2])
        @test ischemicalequal(pt1.Chemical[2], pt1.Chemical[2])
        @test ischemicalequal(chemicalentity(pt1.Chemical[1]), analyzedchemical(pt1.Chemical[1]))
        @test ischemicalequal(inputchemical(pt1.Chemical[1]), pt1.Chemical[1])
        @test ischemicalequal(outputchemical(pt1.Chemical[1]), pt1.Chemical[1])
        @test detectedcharge(pt1.Chemical[1]) == charge(pt1.Chemical[1])
        @test detectedelements(pt1.Chemical[1]) == chemicalelements(pt1.Chemical[1])
        @test detectedisotopes(pt1.Chemical[2]) == isotopomersisotopes(pt1.Chemical[2])
        @test isapprox(retentiontime(pt1.Chemical[1]), 1.5)
        @test isapprox(mz(inputchemical(pt2.Chemical[1])), mz(chemicaltransition(pt2.Chemical[1])[1]))
        @test isapprox(mmi(outputchemical(pt2.Chemical[1])), mmi(detectedchemical(pt2.Chemical[1])))
        @test chemicalentity(analyzedchemical(pt2.Chemical[1])) == chemicalentity(pt2.Chemical[1])
        @test chemicalparent(seriesanalyzedchemical(pt2.Chemical[1])[1]) == chemicalparent(pt2.Chemical[1])
        @test msstage(pt2.Chemical[2]) == 2
    end

    @testset "ChemicalTransition" begin
        @test chemicalname(cp1) == string(chemicalname(inputchemical(cp1)), " -> ", chemicalname(outputchemical(cp1)))
        @test chemicalabbr(cp1) == string(chemicalabbr(inputchemical(cp1)), " -> ", chemicalabbr(outputchemical(cp1)))
        @test chemicalformula(cp1) == chemicalformula(chemicalelements(cp1))
        @test isempty(isotopomersisotopes(cp1))
        @test isapprox(retentiontime(cp1), 7.8)
        @test charge(cp1) == -1
        @test isapprox(mz(cp1), mz(icps[1]))
        @test isapprox(mz(cp1, "[M+H]+"), mz(icps[1], "[M+H]+"))
        @test isapprox(mmi(cp1), mmi(icps[1]))
        @test isapprox(molarmass(cp1), molarmass(icps[1]))
        @test isempty(chemicalsmiles(cp1))
        @test msstage(cp2) == 3
        @test ischemicalequal(cp1, cp1)
        @test ischemicalequal(inputchemical(cp1), first(chemicaltransition(cp1)))
        @test ischemicalequal(outputchemical(cp1), last(chemicaltransition(cp1)))
        @test ischemicalequal(chemicalentity(cp1), chemicalentity(first(chemicaltransition(cp1))))
        @test ischemicalequal(analyzedchemical(cp2), analyzedchemical(cp3))
        @test ischemicalequal(detectedchemical(cp2), detectedchemical(cp3))
        @test ischemicalequal(last(seriesanalyzedchemical(cp2)), detectedchemical(cp2))
        @test ischemicalequal(seriesanalyzedchemical(cp2)[2], icps[1])
        @test ischemicalequal(seriesanalyzedchemical(sp2)[2], ipsi1[1])
        @test ischemicalequal(detectedchemical(sp4), AdductIon(psi2.fa1, ChemicalLoss(Proton())))
        @test ischemicalequal(chemicalparent(cp1), ChemicalTransition(chemicalparent.(chemicaltransition(cp1))...))
        @test detectedisotopes(sp5) == last(seriesanalyzedisotopes(sp5))
        @test detectedisotopes(cp5) == last(seriesanalyzedisotopes(cp5))
        @test detectedelements(sp5) == last(seriesanalyzedelements(sp5))
        @test detectedelements(cp0) == last(seriesanalyzedelements(cp0))
        @test detectedcharge(sp5) == last(seriesanalyzedcharge(sp5))
        @test detectedcharge(cp5) == last(seriesanalyzedcharge(cp5))
        @test detectedisotopes(itit12.Chemical[1]) == last(seriesanalyzedisotopes(itit12.Chemical[1]))
        @test detectedisotopes(itit12.Chemical[11]) == last(seriesanalyzedisotopes(itit12.Chemical[11]))
        @test detectedcharge(itit12.Chemical[11]) == last(seriesanalyzedcharge(itit12.Chemical[11]))
        @test MSC.dictionary_elements(seriesanalyzedelements(cp5)[2]) == MSC.dictionary_elements(MSC.gain_elements(chemicalelements(cp5), chemicalelements(chemicaltransition(cp5)[2])))
        @test MSC.dictionary_elements(seriesanalyzedelements(sp5)[2]) == MSC.dictionary_elements(chemicalelements(isotopomerize(completescheme(ipsi1[1], outputchemical(sp1)), ["[13C]" => 5])))
        @test MSC.dictionary_elements.(chemicalelements.(seriesanalyzedchemical(itit12.Chemical[11]))) == MSC.dictionary_elements.(chemicalelements.(seriesanalyzedchemical(itit12.Chemical[11])))
    end

    @testset "IsotopomerizedSchema" begin
        @test IsotopomerizedSchema(chemicaltransition(itit12.Chemical[2])[2].parent, "-H2OO[13C]") == chemicaltransition(itit12.Chemical[5])[2]
        @test @test_noerror test_show(itit12.Chemical[5])
        @test chemicalabbr(itit12.Chemical[5]) == join(chemicalabbr.(chemicaltransition(itit12.Chemical[5])), " -> ")
        @test chemicalname(itit12.Chemical[5]) == join(chemicalname.(chemicaltransition(itit12.Chemical[5])), " -> ")
        @test chemicalformula(chemicaltransition(itit12.Chemical[5])[2]) == chemicalformula(chemicalelements(chemicaltransition(itit12.Chemical[5])[2]))
        @test MSC.unique_elements(Dict, isotopomersisotopes(chemicaltransition(itit12.Chemical[5])[2])) == MSC.dictionary_elements(filter(!iselement ∘ first, chemicalelements(chemicaltransition(itit12.Chemical[5])[2])))
        @test ischemicalequal(chemicaltransition(itit12.Chemical[1])[2], chemicaltransition(itit12.Chemical[2])[2])
        @test !ischemicalequal(chemicaltransition(itit12.Chemical[1])[2], chemicaltransition(itit12.Chemical[2])[3])
        @test ischemicalequal(chemicaltransition(itit12.Chemical[1])[2], chemicalparent(chemicaltransition(itit12.Chemical[1])[2]))
        @test chemicalparent.(chemicaltransition(itit12.Chemical[5])) == chemicalparent.(chemicaltransition(itit12.Chemical[1]))
        @test elementalscheme(chemicaltransition(itit12.Chemical[1])[2]) == structuralscheme(chemicaltransition(itit12.Chemical[1])[2]) 
        @test isapprox(sum(mmi, chemicaltransition(itit12.Chemical[5])), itit12.MZ3[5]; rtol = 20e-6)
        @test isapprox(sum(molarmass, chemicaltransition(itit12.Chemical[5])), molarmass(detectedchemical(itit12.Chemical[5])); rtol = 20e-6)
    end
    
    @testset "Groupedisotopomerizedschema" begin
        @test @test_noerror test_show(gitit12.Chemical[6])
        @test chemicalabbr(gitit12.Chemical[6]) == join(chemicalabbr.(chemicaltransition(gitit12.Chemical[6])), " -> ")
        @test chemicalname(gitit12.Chemical[6]) == join(chemicalname.(chemicaltransition(gitit12.Chemical[6])), " -> ")
        @test chemicalformula(chemicaltransition(gitit12.Chemical[6])[2]) == chemicalformula(chemicalelements(chemicaltransition(gitit12.Chemical[6])[2]))
        @test MSC.unique_elements(isotopomersisotopes(chemicaltransition(gitit12.Chemical[6])[2])) == MSC.unique_elements(filter(!iselement ∘ first, chemicalelements(chemicaltransition(gitit12.Chemical[6])[2])))
        @test groupedisotopomersabundance(chemicaltransition(gitit12.Chemical[6])[2]) == chemicaltransition(gitit12.Chemical[6])[2].abundance
        @test chemicalparent.(chemicaltransition(gitit12.Chemical[5])) == chemicalparent.(chemicaltransition(gitit12.Chemical[1]))
        @test ischemicalequal(gitit12.Chemical[5], gitit12.Chemical[5])
        @test ischemicalequal(gitit12.Chemical[1], itit12.Chemical[1])
        @test elementalscheme(chemicaltransition(gitit12.Chemical[1])[2]) == structuralscheme(chemicaltransition(gitit12.Chemical[1])[2]) 
        @test isapprox(sum(mmi, chemicaltransition(gitit12.Chemical[6])), gitit12.MZ3[6]; rtol = 20e-6)
        @test isapprox(sum(molarmass, chemicaltransition(gitit12.Chemical[6])), molarmass(detectedchemical(gitit12.Chemical[6])); rtol = 20e-6)
    end
end
