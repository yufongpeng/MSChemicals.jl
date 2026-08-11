@info "Running tests for Spectrum and CoelutingIsobars"

@testset "MS analyzer" begin 
    @test all((q.window isa W for (q, W) in zip(qa, [MSC.SuperGaussianWindow, MSC.GaussianWindow, MSC.CosineWindow, MSC.TukeyWindow, MSC.RectWindow])))
    @test all((q.window isa W for (q, W) in zip(qita, [MSC.GaussianWindow, MSC.SuperGaussianWindow])))
    @test all((q.window isa W for (q, W) in zip(lita, [MSC.GaussianWindow, MSC.CosineWindow, MSC.TukeyWindow, MSC.RectWindow])))
    for q in qa 
        dw = MSC.discrete_window(q, [500.0], 0.0001, 1, 0.0001)
        @test isapprox(0.0001 * (findlast(>=(0.5), first(dw)) - findfirst(>=(0.5), first(dw))), MSC.fwhm_mz(q, 500.0); rtol = 1e-3)
        @test isapprox(q.window(500.0, MSC.window_parameter(q, 500.0)...), 1.0)
        @test isapprox(q.window(500.35, MSC.window_parameter(q, 500.0)...), 0.5) || isapprox(q.window(500.7, MSC.window_parameter(q, 500.0)...), 0.5) || isapprox(q.window(500.35, MSC.window_parameter(q, 500.0)...), 0.0)
    end
    for q in msa 
        dw = MSC.discrete_window(q, [500.0], 0.0001, 1, 0.0001)
        @test isapprox(0.0001 * (findlast(>=(0.5), first(dw)) - findfirst(>=(0.5), first(dw))), MSC.fwhm_mz(q, 500.0); rtol = 1e-3)
        @test isapprox(q.window(500.0, MSC.window_parameter(q, 500.0)...), 1.0)
        @test isapprox(q.window(500.5, MSC.window_parameter(q, 500.0)...), 0.5) || isapprox(q.window(500.5, MSC.window_parameter(q, 500.0)...), 0.0)
    end
    for (ms, nm) in [Quadrupole() => "Quadrupole", QIT() => "Quadrupole Ion Trap", LIT() => "Linear Ion Trap", TOF() => "TOF", Orbitrap() => "Orbitrap", FTICR() => "FTICR"]
        @test MSC.msanalyzer_name(ms) == nm
    end
    for (ms, nm) in zip(msa, ["MS-Analyzer with Uniform Window (Gaussian-Tailed)", "MS-Analyzer with Power Cosine Window", "MS-Analyzer with Tukey Window (Fixed Taper)", "MS-Analyzer with Rectangular Window", "MS-Analyzer with Cosine Window"])
        @test MSC.msanalyzer_name(ms) == nm
    end
    @test @test_noerror plot_window(MSC.GaussianWindow())
    @test @test_noerror plot_window!(MSC.TukeyWindow(0.2))
    @test @test_noerror plot_resolving_power((0, 1000), FTICR())
    @test @test_noerror plot_resolving_power!((0, 1000), Orbitrap())
    @test @test_noerror plot_resolving_power(QuadrupoleIonTrap((0, 1000)))
    @test @test_noerror plot_resolving_power!(LinearIonTrap((0, 1000)))
end

@testset "Spectrum" begin 
    @test @test_noerror test_show(spec1)
    @test isapprox(pt2.Abundance2[1], 20000)
    @test @test_noerror plot_spectrum((948.1, 951.5), spec2)
    @test @test_noerror plot_spectrum!((948.1, 951.5), spec2; deconvolution = true)
    @test @test_noerror plot_spectrum(spec2)
    @test @test_noerror plot_spectrum!(spec2)
    @test @test_noerror plot_spectrum!(spec2.table)
end

@testset "CoelutingIsobars" begin 
    @test @test_noerror test_show(ci1)
    @test retentiontime(ci1.target.Chemical[4]) - retentiontime(ci1.tables[4].Chemical[1]) <= 0.15
    @test isapprox(sum(ci1.tables[4].Abundance1[2:3]), ib1.var"Abundance1(%)"[4] * ci1.target.Abundance1[4] / 100)
    @test all(x -> abs(x) < 0.1, ib2.var"ΔMZ1")
    @test all(>(100 * 1e-4), ib2.var"Abundance2(%)")
end
