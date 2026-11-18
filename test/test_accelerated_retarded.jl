function accelerated_remainder_run(N)
    q, z0, A, ω = 0.7, 0.0, 0.08, 1.6
    t0, tf, L = 0.70, 0.78, 2.0
    zfun(t) = z0 + A * sin(ω * t)
    vfun(t) = A * ω * cos(ω * t)
    x, h = centered_grid(L, N)
    params = forced_params(x, h, q, z0, A, ω)

    initial = [ref_remainder(xi, t0, q, zfun, vfun) for xi in x]
    Q0 = hcat(getproperty.(initial, :Π),
              getproperty.(initial, :Ψ),
              getproperty.(initial, :Φ))
    dt = 0.2h
    sol = solve(ODEProblem(ODE.field_rhs_forced!, Q0, (t0, tf), params), SSPRK54();
                adaptive=false, dt=dt, saveat=[tf])
    @test sol.retcode == ReturnCode.Success
    Q = sol.u[end]

    exact = [ref_remainder(xi, tf, q, zfun, vfun) for xi in x]
    mask = abs.(x) .< 0.45
    Πe = getproperty.(exact, :Π)
    Ψe = getproperty.(exact, :Ψ)
    Φe = getproperty.(exact, :Φ)
    field_error = sqrt(h * sum((Q[mask, 1] .- Πe[mask]).^2 .+
                               (Q[mask, 2] .- Ψe[mask]).^2 .+
                               (Q[mask, 3] .- Φe[mask]).^2))
    _, Πp, Ψp = ODE.regular_fields_at_particle(Q, zfun(tf), params)
    force_error = abs(Ψp + vfun(tf) * Πp)
    return field_error, force_error
end

@testset "accelerated retarded-minus-puncture convergence" begin
    errors = accelerated_remainder_run.((100, 200, 400))
    field_errors = first.(errors)
    force_errors = last.(errors)
    @test field_errors[2] < field_errors[1]
    @test field_errors[3] < field_errors[2]
    @test force_errors[2] < force_errors[1]
    @test force_errors[3] < force_errors[2]

    q, A, ω, t = 0.7, 0.08, 1.6, 0.73
    zfun(s) = A * sin(ω * s)
    vfun(s) = A * ω * cos(ω * s)
    for side in (-1.0, 1.0), ϵ in (1e-4, 1e-6)
        rem = ref_remainder(zfun(t) + side * ϵ, t, q, zfun, vfun)
        @test abs(rem.Φ) < 2e-4
        @test abs(rem.Π) < 2e-4
        @test abs(rem.Ψ) < 2e-4
    end
end
