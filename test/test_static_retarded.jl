@testset "exact stationary retarded solution" begin
    q, m0, t = 0.8, 3.0, 1.25
    zfun(s) = 0.0
    vfun(s) = 0.0

    for x in (-1.4, -0.7, 0.3, 1.4)
        exact = ref_retarded_fields(x, t, q, zfun, vfun)
        params = (z0=0.0, A=0.0, ω=1.0, q=q)
        @test RP.Φ_ret(x, params, t) ≈ exact.Φ
        @test RP.Π_ret(x, params, t) ≈ exact.Π
        @test RP.Ψ_ret(x, params, t) ≈ exact.Ψ
        if abs(x) < t
            @test exact.Φ ≈ -q * (t - abs(x)) / 2
            @test exact.Π ≈ -q / 2
            @test exact.Ψ ≈ q * sign(x) / 2
            rem = ref_remainder(x, t, q, zfun, vfun)
            @test rem.Φ ≈ 0 atol=1e-13
            @test rem.Π ≈ 0 atol=1e-13
            @test rem.Ψ ≈ 0 atol=1e-13
        else
            @test exact == (Φ=0.0, Π=0.0, Ψ=0.0)
        end
    end

    m = m0 - q^2 * t / 2
    Efield = q^2 * t / 2
    @test m + Efield ≈ m0
end
