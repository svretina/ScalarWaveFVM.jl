ref_gamma(v) = inv(sqrt(1 - v^2))

function ref_puncture(x, z, v, q, τ)
    γ = ref_gamma(v)
    r = x - z
    n = sign(r)
    return (Φ=q * (γ * (abs(r) + v * r) - τ) / 2,
            Π=-q * γ * (1 + v * n) / 2,
            Ψ=q * γ * (n + v) / 2)
end

function ref_retarded_time(x, t, zfun)
    t <= 0 && return nothing
    f(tr) = t - tr - abs(x - zfun(tr))
    f(0.0) <= 0 && return nothing
    return find_zero(f, (0.0, t), Bisection())
end

function ref_retarded_fields(x, t, q, zfun, vfun)
    tr = ref_retarded_time(x, t, zfun)
    tr === nothing && return (Φ=0.0, Π=0.0, Ψ=0.0)
    τr, _ = quadgk(s -> sqrt(1 - vfun(s)^2), 0.0, tr)
    zr, vr = zfun(tr), vfun(tr)
    n = sign(x - zr)
    invγr = sqrt(1 - vr^2)
    return (Φ=-q * τr / 2,
            Π=-q * invγr / (2 * (1 - n * vr)),
            Ψ=q * n * invγr / (2 * (1 - n * vr)))
end

function ref_remainder(x, t, q, zfun, vfun)
    τ, _ = quadgk(s -> sqrt(1 - vfun(s)^2), 0.0, t)
    ret = ref_retarded_fields(x, t, q, zfun, vfun)
    puncture = ref_puncture(x, zfun(t), vfun(t), q, τ)
    return (Φ=ret.Φ - puncture.Φ,
            Π=ret.Π - puncture.Π,
            Ψ=ret.Ψ - puncture.Ψ)
end

function centered_grid(L, N)
    h = 2L / N
    return collect(range(-L + h / 2; step=h, length=N)), h
end

function forced_params(x, h, q, z0, A, ω)
    return (equation=EQ.LinearScalarWaveEquation1D(1.0),
            numerical_flux=NF.flux_godunov,
            slope_limiter=LIM.MC,
            h=h, N=length(x), x=x, q=q, z0=z0, A=A, ω=ω)
end
