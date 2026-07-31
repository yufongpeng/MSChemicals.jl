@info "Defining MS analyzer"

qa = [
    Quadrupole(; fwhm = 1.4),
    Quadrupole(; fwhm = 0.7), 
    Quadrupole(500.0; offset = 0.1, taperproportion = 1),
    Quadrupole(500.0; taperproportion = 0.5),
    Quadrupole(500.0; taperproportion = 0)
]

qita = [
    QIT(),
    QIT(500.0), 
]

lita = [
    LIT(),
    LIT(500.0; offset = 0.1, taperproportion = 1),
    LIT(500.0; taperproportion = 0.5),
    LIT(500.0; taperproportion = 0)
]

msa = [
    MSAnalyzer(MSC.GaussianTailedUniformWindow(0.5), 500.0, 0.1, x -> 1.0, 0.1, 0.1),
    MSAnalyzer(MSC.PowerCosineWindow(0.5), 500.0, 0.1, x -> 1.0, 0.1, 0.1),
    MSAnalyzer(MSC.FixedTaperTukeyWindow(0.2), 500.0, 0.1, x -> 1.0, 0.1, 0.1),
    MSAnalyzer(MSC.RectWindow(), 500.0, 0.1, x -> 1.0, 0.1, 0.1),
    MSAnalyzer(MSC.CosineWindow(), 500.0, 0.1, x -> 1.0, 0.1, 0.1)
]

@info "Defining chemical and product table"

GQ1 = Chemical("GQ1b 18:1;O2/18:0", "C106H182N6O55")
GD3 = Chemical("GD3 18:1;O2/18:0", "C70H125N3O29")
nana = Chemical("N-acetylneuraminic acid", "C11H19NO9")
galnac = Chemical("GalNAc", "C8H15NO6")
f366 = Chemical("GalGalNAc", "C6H12O6C8H13NO5")
f657 = Chemical("NeuAcGalGalNAc", "C11H19NO9C6H10O5C8H13NO5")
f948 = Chemical("NeuAc2GalGalNAc", "C22H36N2O17C6H10O5C8H13NO5")
pt = Table(; 
    Chemical = [
        AdductIon(GQ1, "[M+2H]2+"), 
        AdductIon(GQ1, "[M+H+Na]2+"), 
        AdductIon(GQ1, "[M+2H-H2O]2+"), 
        AdductIon(GD3, "[M+H]+"), 
        AdductIon(f948, "[M+H-H2O]+"), 
        AdductIon(f657, "[M+H-H2O]+"), 
        AdductIon(f366, "[M+H-H2O]+"), 
        AdductIon(galnac, "[M+H-H2O]+"), 
        AdductIon(nana, "[M+H-H2O]+"), 
        AdductIon(nana, "[M+H-2H2O]+")], 
    Product = [
        [
            AdductIon(f948, "[M+H-H2O]+"), 
            AdductIon(f657, "[M+H-H2O]+"), 
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-2H2O]+"),
            AdductIon(GQ1, "[M+2H]2+"), 
            AdductIon(GQ1, "[M+2H-H2O]2+"), 
            AdductIon(GD3, "[M+H]+")
        ],
        [
            AdductIon(f948, "[M+H-H2O]+"), 
            AdductIon(f657, "[M+H-H2O]+"), 
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+"),
            AdductIon(GQ1, "[M+H+Na]2+"), 
            AdductIon(GD3, "[M+H]+")
        ],
        [
            AdductIon(f948, "[M+H-H2O]+"), 
            AdductIon(f657, "[M+H-H2O]+"), 
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+"),
            AdductIon(GQ1, "[M+2H-H2O]2+"), 
            AdductIon(GD3, "[M+H]+")
        ],
        [
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+"),
            AdductIon(GD3, "[M+H]+")
        ],
        [
            AdductIon(f948, "[M+H-H2O]+"), 
            AdductIon(f657, "[M+H-H2O]+"), 
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+")
        ],
        [
            AdductIon(f657, "[M+H-H2O]+"), 
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+")
        ],
        [
            AdductIon(f366, "[M+H-H2O]+"), 
            AdductIon(galnac, "[M+H-H2O]+"), 
        ],
        [
            AdductIon(galnac, "[M+H-H2O]+")
        ],
        [
            AdductIon(nana, "[M+H-H2O]+"),
            AdductIon(nana, "[M+H-2H2O]+")
        ],
        [
            AdductIon(nana, "[M+H-2H2O]+")
        ],
    ], 
    Proportion = [
        [0.3, 0.1, 0.1, 0.05, 0.1, 0.2, 0.3, 0.05, 0.1], 
        [0.35, 0.12, 0.12, 0.07, 0.1, 0.2, 0.2, 0.2], 
        [0.35, 0.12, 0.2, 0.07, 0.12, 0.3, 0.2, 0.2], 
        [0.1, 0.4, 0.5], 
        [0.05, 0.15, 0.15, 0.3, 0.1, 0.3], 
        [0.05, 0.1, 0.4, 0.15, 0.4], 
        [0.1, 0.7], 
        [0.8],
        [0.2, 0.8], 
        [0.8]
    ])
IGQ1b = AdductIon(GQ1, "[M+2H]2+")
IGD3 = AdductIon(GD3, "[M+H]+")
ms0 = Table(; Chemical = [GQ1, GD3], Adduct = [["[M+2H]2+", "[M+H+Na]2+", "[M+2H-H2O]2+"], ["[M+H]+"]], Proportion = [[1, 0.1, 0.3], [1]], Abundance = [100000, 10000])
ms1 = Table(; Chemical = repeat([GQ1, GD3], 1000))

@info "Running HRMS and MS/MS"

Ionization(repeat([GQ1, GD3], 1000); abundance = 10, proportion = [100, 10], adduct = ["[M+2H]2+", "[M+2H-H2O]2+"])
Ionization(ms1; adduct = ["[M+2H]2+", "[M+2H-H2O]2+"])

spec1 = @p ms0 |> Ionization(; abtype = :list) |> MSScan
spec2 = @p spec1 |> Isolation(Quadrupole(1210.588728; fwhm = 1.3, offset = 0.3)) |> Fragmentation(pt; threading = false) |> MSScan
spec3 = @p spec2 |> Isolation(Quadrupole(948.3303; fwhm = 1.3, offset = 0.3)) |> Fragmentation(pt; threading = true) |> MSScan
spec4 = @p spec3 |> AllIons((1000, 2000)) |> MSScan
spec5 = @p spec3 |> Fragmentation(pt) |> MSScan
spec6 = @p spec5 |> AllIons() |> MSScan
spec7 = @p spec6 |> MSScan(Quadrupole())

@info "Running SIM"

transitiontable = Table(;
    Transition = collect(1:5),
    MS1 = [
            Quadrupole(mz(IGQ1b); fwhm = 0.7, digits = 1), 
            Quadrupole(mz(IGQ1b); fwhm = 0.7, digits = 1), 
            Quadrupole(mz(IGQ1b); fwhm = 0.7, digits = 1),
            Quadrupole(mz(IGD3); fwhm = 0.7, digits = 1),
            Quadrupole(1000.0; fwhm = 0.7, digits = 1)
    ],
    Product = [pt, pt, pt, pt, pt],
    MS2 = [
        Quadrupole(mz(AdductIon(nana, "[M+H-2H2O]+")); fwhm = 0.7, digits = 1),
        Quadrupole(mz(AdductIon(f948, "[M+H-H2O]+")); fwhm = 0.7, digits = 1), 
        Quadrupole(mz(IGD3); fwhm = 0.7, digits = 1),
        Quadrupole(mz(AdductIon(nana, "[M+H-2H2O]+")); fwhm = 0.7, digits = 1),
        Quadrupole(200.0; fwhm = 0.7, digits = 1)
    ]
)
sim = SelectedIonMonitor(transitiontable, MSScan(vcat(Isotopologues(IGQ1b; abtype = :input, abundance = 100000), Isotopologues(IGD3; abtype = :input, abundance = 10000))))
pt2 = peak_table(sim)

@info "Running CoelutingIsobars"

exp1 = Table(;
    Chemical = AdductIon.([
        Chemical("Fructose", "C6H12O6"; retentiontime = 1.25),
        Chemical("Glucose", "C6H12O6"; retentiontime = 1.51), 
        Chemical("Galactose", "C6H12O6"; retentiontime = 2.005), 
        Chemical("Mannose", "C6H12O6"; retentiontime = 1.97), 
        Chemical("Glucosamine", "C6H13O5N"; retentiontime = 2.12),
        Chemical("Galactosamine", "C6H13O5N"; retentiontime = 2.61),
        Chemical("Mannosamine", "C6H13O5N"; retentiontime = 2.51),
    ], "[M+H]+"),
    Abundance1 = [100, 10000, 5000, 4000, 2000, 1000, 700], 
)

ci1 = CoelutingIsobars(
        [retentiontime => acrit(0.15)], 
        [
            TOF([]; resolution = 20000) => rcrit(1e-5)
        ], 
        exp1
    )
ib1 = isobar_table(ci1; threshold = rcrit(1e-5))

FA180 = AdductIon(Chemical("FA 18:0", "C18H36O2"), "[M-H]-")
FA181 = AdductIon(Chemical("FA 18:1", "C18H34O2"), "[M-H]-")
FA200 = AdductIon(Chemical("FA 20:0", "C20H40O2"), "[M-H]-")
exp2 = Table(;
    Chemical = ChemicalSeries.([
        AdductIon(Chemical("PC O-18:1/18:0", "C44H88NO7P"; retentiontime = 7.8), "[M+OFo]-") => FA180,
        AdductIon(Chemical("PC O-18:0/18:1", "C44H88NO7P"; retentiontime = 7.85), "[M+OFo]-") => FA181,
        AdductIon(Chemical("PS 18:0/20:0", "C44H86NO10P"; retentiontime = 7.7), "[M-H]-") => FA180,
        AdductIon(Chemical("PS 18:0/20:0", "C44H86NO10P"; retentiontime = 7.7), "[M-H]-") => FA200,
        AdductIon(Chemical("PC P-18:0/18:0", "C44H88NO7P"; retentiontime = 7.94), "[M+OFo]-") => FA180,
        AdductIon(Chemical("PC P-18:0/18:1", "C44H86NO7P"; retentiontime = 7.68), "[M+OFo]-") => FA181,
    ]),
    Abundance1 = [2000, 1500, 10000, 5000, 4000, 3000], 
)

ci2 = CoelutingIsobars(
        [retentiontime => acrit(0.15)], 
        [
            Quadrupole([]; accuracy = 0.1, fwhm = 1.3, digits = 1) => rcrit(1e-3),
            Quadrupole([]; accuracy = 0.1, fwhm = 1.3, digits = 1) => rcrit(1e-3)
        ], 
        exp2
    )
ib2 = isobar_table(ci2)