@info "Running tests for Utils"

@testset "Utils" begin 
    @testset "Combinatorics" begin 
        @test MSC.safe_factorial(Val(false), 1, 1) == 1
        @test MSC.safe_factorial(1, 1) == 1
        @test isapprox(MSC.safe_factorial(50, 50), factorial(50, 50))
        @test isapprox(MSC.safe_factorial(Val(false), 10, 3), factorial(10, 3))
        @test isapprox(MSC.safe_factorial(Val(false), 50, 3), MSC.safe_factorial(Val(true), 50, 3))
        @test isapprox(MSC.safe_multinomial(Val(false), [50, 50]), MSC.safe_multinomial(Val(false), 50, 50))
        @test isapprox(MSC.safe_multinomial(Val(true), [50, 50]), MSC.safe_multinomial(Val(true), 50, 50))
        @test MSC.safe_multinomial(Val(false), 50) == MSC.safe_multinomial(Val(true), 50)
        @test MSC.safe_multinomial(50) == 1
        @test isapprox(MSC.safe_multinomial(50, 50), MSC.safe_multinomial([50, 50]))
        @test isapprox(MSC.safe_multinomial(Val(false), [10, 10]), MSC.safe_multinomial(Val(false), 10, 10))
    end
    @testset "New elements" begin 
        m = [23.98504168, 24.985836966, 25.982592972]
        a = [0.78965, 0.1001, 0.11025]
        set_element!("Mg", m , a)
        @test isapprox(mmi("Mg"), m[1])
        @test isapprox(molarmass("Mg"), m'a)
        @test iselement("Mg")
        @test isisotope("[25Mg]")
        @test MSC.ismajor("[24Mg]")
        @test MSC.isminor("[26Mg]")
    end
    @testset "Criteria and Interval" begin 
        @test @test_noerror test_show(ri"[0, 10]")
        @test @test_noerror test_show(crit(ri"()", ri"[0, 10]"))
        @test @test_noerror test_show(acrit(ri"[0, 10]"))
        @test @test_noerror test_show(rcrit(ri"[0, 10]"))
        @test !qualified_peak1(85, 100, ct1)
        @test qualified_peak1(85, 100, ct2)
        @test !qualified_peak2(8, 40, ct1)
        @test qualified_peak2(8, 40, ct2)
        @test union(ri"[1,3)", ri"[2,4)") == ri"[1,4)"
        @test 2 in ri"(-Inf, 5]"
        @test ri"(-5, 5]" * (-10) / 2 + 10 - 5 == ri"[-20, 30)"
        @test ri"(-∞, ∞)" * (-10) / 2 + 10 - 5 == ri"(-∞, ∞)"
        @test ri"(-5, 5]" * Inf64 == ri"(-5, 5]" / 0
        @test ri"(-5, 5]" + Inf64 == ri"(-5, 5]" - Inf64
    end
    @testset "Measure" begin 
        @test isapprox(relative_error(10, 100) * 100, percentage_error(10, 100))
        @test isapprox(relative_error_mean(10, 100) * 100, percentage_error_mean(10, 100))
        @test isapprox(relative_error(0.01, 100) * 1e6, ppm_error(0.01, 100))
        @test isapprox(relative_error_mean(0.01, 100) * 1e6, ppm_error_mean(0.01, 100))
        @test MSC.measure_name(charge, value_error) == MSC.measure_name("Z", value_error)
        @test MSC.measure_name(molarmass, relative_error) == MSC.measure_name("M", relative_error_mean)
        @test MSC.measure_name(mz, percentage_error) == MSC.measure_name("MZ", percentage_error_mean)
        @test MSC.measure_name(mmi, ppm_error) == MSC.measure_name("Mmi", ppm_error_mean)
        @test MSC.measure_name(sum, relative_error) == "Δsum/sum"
    end
    @testset "Others" begin 
        @test MSC.lastcolnum(["MZ1", "Abundance2", "MZ3", "MZ2", "Abundance1", "Abundance3"], "MZ") == :MZ3
        @test MSC.ithcolnum(["MZ1", "Abundance2", "MZ3", "MZ2", "Abundance1", "Abundance3"], "MZ", 1) == :MZ1
        @test MSC.allcolnum(["MZ1", "Abundance2", "MZ3", "MZ2", "Abundance1", "Abundance3"], "Abundance") == [:Abundance1, :Abundance2, :Abundance3]
        @test MSC.normalize_abundance([1, 2, 3, 4], 1, MSC.abtyped(:raw)) == [1, 2, 3, 4]
        @test MSC.normalize_abundance([1, 2, 3, 4], 1, MSC.abtyped(:total)) == [1, 2, 3, 4]
        @test MSC.normalize_abundance([1, 2, 3, 4], 1, MSC.abtyped(:input)) == [1, 2, 3, 4] ./ 1
        @test all(isapprox.(MSC.normalize_abundance([1, 2, 3, 4], 1, MSC.abtyped(MSC.deabtyped(MSC.Max()))), [1, 2, 3, 4] ./ 4))
        @test all(isapprox.(MSC.normalize_abundance([1, 2, 3, 4], 1, MSC.abtyped(MSC.deabtyped(:list))), [1, 2, 3, 4] ./ 10))
        @test match_chemical(AdductIon.([Chemical("Fructose", "C6H12O6"), Chemical("Glucose", "C6H12O6"), Chemical("Galactose", "C6H12O6")], "[M+H]+"), exp1).LibID == [1, 2, 3]
    end
end