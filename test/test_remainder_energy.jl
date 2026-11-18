@testset "continuum remainder energy identity" begin
    q, v, a = 0.8, -0.35, 0.27
    γ = ref_gamma(v)
    x = collect(range(-1.0, 1.0; length=2001))
    h = x[2] - x[1]
    Πr = @. 0.2 + 0.1 * x
    Ψr = @. -0.15 + 0.05 * x
    n = sign.(x)
    SΠ = @. q * a * γ^3 * (v + n) / 2
    SΨ = @. -q * a * γ^3 * (1 + v * n) / 2

    source_work = EN.newton_cotes(Πr .* SΠ .+ Ψr .* SΨ, h; order=1)
    boundary_term = Πr[end] * Ψr[end] - Πr[1] * Ψr[1]
    dΠ = 0.05 .+ SΠ # ∂xΨr + SΠ
    dΨ = 0.1 .+ SΨ # ∂xΠr + SΨ
    direct_rate = EN.newton_cotes(Πr .* dΠ .+ Ψr .* dΨ, h; order=1)
    @test direct_rate ≈ boundary_term + source_work atol=5e-7
end
