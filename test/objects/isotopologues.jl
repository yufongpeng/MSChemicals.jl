@info "Running Isotopologues and TandemIsotopologues"

p1 = MSC.maximal_proportion(Val(false), Dict("[33S]" => 6, "[34S]" => 6), Dict("S" => 1))
p2 = MSC.maximal_proportion(Val(false), Dict("[33S]" => 6, "[34S]" => 6), Dict("S" => 11))

it1 = Isotopologues(icglc[1]; abundance = 1e5, threshold = crit(1e1, 1e-2))
it2 = Isotopologues(ioncore(icglc[1]); abundance = 1e5, abtype = :total, threshold = crit(1e1, 1e-2))
it3 = Isotopologues(ipsi2[1]; abtype = :total, threshold = crit(1e-3, 1e-3))
it4 = Isotopologues("C6H12O6"; abundance = 1e5, abtype = :total, threshold = crit(1e1, 1e-2))
it5 = Isotopologues(Table(; Chemical = repeat([icglc[1]], 1000), Abundance1 = repeat([1e5], 1000)); threshold = crit(1e1, 1e-2), threading = true)
it6 = Isotopologues(repeat([icglc[1]], 1000); abundance = 1e5, threshold = crit(1e1, 1e-2), threading = false)

git3 = group_isotopologues(it3)

# Tandemisotopologues / Isotopologues MS/MS
itit1 = TandemIsotopologues(cp1; abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.3)
itit2 = TandemIsotopologues(chemicalformula(icps[1]) => chemicalformula(last(chemicaltransition(cp1))); abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.3, chemicalparser = ChemicalTransitionParser())
itit3 = TandemIsotopologues(ChemicalSeries(AdductIon(cps, "[M-H]-"), lossserine); abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.5)
itit4 = TandemIsotopologues(string("[", chemicalformula(cps), "-H]-") => "-C3H5NO2"; abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.5)
itit5 = TandemIsotopologues(ChemicalSeries(ips[1], AdductIon(ioncore(ips[1]).fa1, "[M-H]-")); abtype = :total, threshold = crit(1e-8, 1e-8), transmission = 0.3)
itit6 = TandemIsotopologues(ChemicalSeries(ipsi1[3], LossSerine(), SN1Acyl()); abtype = :total, abundance = 1e5, threshold = crit(1e1, 1e-2), transmission = 0.5, precise = true)

itit7 = Isotopologues("[S8]2-" => "[S7]-"; abtype = :total)
itit8 = Isotopologues("[C2H5[13C]OO]-" => "[C2H5]-"; abtype = :total)
itit9 = Isotopologues("[C2H5[13C]OO]-" => "-[13C]O2"; abtype = :total)
itit10 = Isotopologues("[C2H5[12C]OO]-" => "-[12C]O2"; abtype = :total)

itit11 = Isotopologues("[C2H5[12C]OO]-" => "[-H2O-CO]" => "+H2O"; abtype = :input)
itit12 = TandemIsotopologues("[C2H5[12C]OO]-" => "[-H2O-CO]" => "+H2O"; abtype = :input)
itit13 = TandemIsotopologues("[C2H5[12C]OO]-" => "-[CO]+" => "+CO2"; abtype = :input)
itit14 = TandemIsotopologues("[C2H5[12C]OO]-" => "-CO" => "+[CO2]-"; abtype = :input)
gitit12 = group_isotopologues(itit12)
gitit13 = group_isotopologues(itit13)
gitit14 = group_isotopologues(itit14)

itit15 = TandemIsotopologues(repeat([icps[3]], 1000); product = [lossserine], abundance = [1e5], threshold = crit(1e1, 1e-2), proportion = [0.5], transmission = 0.3, threading = true)
itit16 = TandemIsotopologues(Table(; Chemical = repeat([icps[3]], 1000), Product = repeat([[lossserine]], 1000), Abundance1 = repeat([1e5], 1000), Proportion = repeat([[0.5]], 1000), Transmission = repeat([0.3], 1000)); threshold = crit(1e1, 1e-2), threading = true)
itit17 = TandemIsotopologues("[C2H5[12C]OO]-" => "-CO" => "-H2O"; abtype = :input)
itit18 = TandemIsotopologues("[C2H5[12C]OO]-" => "-CO" => "OH"; abtype = :input)
itit19 = TandemIsotopologues("[C200H500[12C]O200]-" => "-C100O100" => "-O50H50"; abtype = :input)
itit20 = TandemIsotopologues("[C200H500[12C]O200]-" => "-C100O100" => "[C100H450[12C]O50]-"; abtype = :input)
itit21 = TandemIsotopologues("[C200H500[12C]O200]-" => "[C100H500[12C]O100]-" => "-O50H50"; abtype = :input)
itit22 = TandemIsotopologues("[C200H500[12C]O200]-" => "[C100H500[12C]O100]-" => "[C100H450[12C]O50]-"; abtype = :input)

d0 = MSC.dictionary_elements(chemicalelements(ipsi2[1]))
d1 = MSC.dictionary_elements(chemicalelements(inputchemical(itit5.Chemical[14])))
d2 = MSC.dictionary_elements(chemicalelements(outputchemical(itit5.Chemical[14])))

d3, d4, d5 = last.(MSC.serieschemicaldata(itit6.Chemical[12]))

@info "Running isotpologues function on large chemical"
Isotopologues("C494H776O148N136S4"; abtype = :total)
ti = time()
itl = Isotopologues("C494H776O148N136S4"; abtype = :total)
te = time()
if te - ti > 1
    @info "`Isotopologues` takes too much time for chemical ≈ 10 kDa. Skip tests for larger chemicals."
else
    ti = time()
    itl = Isotopologues("C1482H2328O444N408S12"; abtype = :total)
    te = time()
    if te - ti > 2
        @info "`Isotopologues` takes too much time for chemical ≈ 30 kDa. Skip tests for larger chemicals."
    else
        ti = time()
        itl = Isotopologues("C2964H4656O888N816S24"; abtype = :total)
        te = time()
        if te - ti > 5
            @info "`Isotopologues` takes too much time for chemical ≈ 50 kDa. Skip tests for larger chemicals."
        else
            itl = Isotopologues("C4940H7760O1480N1360S40"; abtype = :total)
        end
    end
end

TandemIsotopologues("C297H388O74N68S2" => "C124H194O37N34S"; abtype = :total)
Fragmentation(Table(; Chemical = ["[C297H388O74N68S2]+"], Product = [["[C124H194O37N34S]+"]]), Isotopologues("[C297H388O74N68S2]+"))
TandemIsotopologues("C297H388O74N68S2" => "C124H194O37N34S" => "C62H87O19N17S"; abtype = :total)
TandemIsotopologues("C297H388O74N68S2"; product = ["C124H194O37N34S"], abtype = :total)
Isotopologues("C297H388O74N68S2" => "C124H194O37N34S"; abtype = :total)
ti = time()
itll1 = TandemIsotopologues("C297H388O74N68S2" => "C124H194O37N34S"; abtype = :total)
itll2 = TandemIsotopologues("C297H388O74N68S2"; product = ["C124H194O37N34S"], abtype = :total)
itll3 = Isotopologues("C297H388O74N68S2" => "C124H194O37N34S"; abtype = :total)
te = time()
if te - ti > 1
    @info "`Isotopologues` and `TandemIsotopologues` takes too much time for transitions ≈ 5 kDa -> 3 kDa. Skip tests for larger transitions."
else
    Fragmentation(Table(; Chemical = ["[C494H776O148N136S4]+"], Product = [["[C247H388O74N68S2]+"]]), Isotopologues("[C494H776O148N136S4]+"))
    ti = time()
    itll1 = TandemIsotopologues("C494H776O148N136S4" => "C247H388O74N68S2"; abtype = :total)
    itll2 = TandemIsotopologues("C494H776O148N136S4"; product = ["C247H388O74N68S2"], abtype = :total)
    itll3 = Isotopologues("C494H776O148N136S4" => "C247H388O74N68S2"; abtype = :total)
    te = time()
    if te - ti > 2
        @info "`Isotopologues` takes too much time for transitions ≈ 10 kDa -> 5 kDa. Skip tests for larger chemicals and transitions."
    else
        Fragmentation(Table(; Chemical = ["[C988H1552O296N272S8]+"], Product = [["[C494H776O148N136S4]+"]]), Isotopologues("[C988H1552O296N272S8]+"))
        ti = time()
        itll1 = TandemIsotopologues("C988H1552O296N272S8" => "C494H776O148N136S4"; abtype = :max)
        itll2 = TandemIsotopologues("C988H1552O296N272S8"; product = ["C494H776O148N136S4"], abtype = :max)
        itll3 = Isotopologues("C988H1552O296N272S8" => "C494H776O148N136S4"; abtype = :max)
        te = time()
        if te - ti > 5
            @info "`Isotopologues` takes too much time for transitions ≈ 20 kDa -> 10 kDa. Skip tests for larger chemicals and transitions."
        else
            itll1 = TandemIsotopologues("C1482H2328O444N408S12" => "C741H1164O222N204S6"; abtype = :max)
            itll2 = TandemIsotopologues("C1482H2328O444N408S12"; product = ["C741H1164O222N204S6"], abtype = :max)
            itll3 = Isotopologues("C1482H2328O444N408S12" => "C741H1164O222N204S6"; abtype = :max)
        end
    end
end