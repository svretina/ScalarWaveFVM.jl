module RetardedPotentials

using NonlinearSolve
using Roots
using QuadGK: quadgk
using ForwardDiff: derivative
using ..ScalarField

export Φ_ret, Π_ret, Ψ_ret, Φ_ret_exact, Π_ret_exact, Ψ_ret_exact

@inline function z(t, p)
    return p.z0 + p.A * sin(p.ω * t)
end

@inline function v(t, p)
    return p.ω * p.A * cos(p.ω * t)
end

@inline function a(t, p)
    return -p.ω^2 * p.A * sin(p.ω * t)
end

function Φ_ret(x, p, t)
    # calculate retarded time
    f(tr, p) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem(f, zero(t), p)
    tᵣ = solve(prob; reltol=1e-18, abstol=1e-18).u
    if tᵣ < 0
        tᵣ = zero(tᵣ)
    end

    z_tr = z(tᵣ, p)
    v_tr = v(tᵣ, p)
    return ScalarField.Φs(x, z_tr, v_tr, p.q, tᵣ)
end

function Φ_ret_exact(x, p, t)
    f(tr) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem((tr, _) -> f(tr), zero(t))
    tr = solve(prob; reltol=1e-18, abstol=1e-18).u
    tr = max(zero(tr), tr)

    integrand(τ) = sqrt(1 - v(τ, p)^2)
    I, _ = quadgk(integrand, zero(tr), tr; rtol=1e-18, atol=1e-18)

    return -p.q * I / 2
end

function Π_ret(x, p, t)
    # calculate retarded time
    f(tr, p) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem(f, zero(t), p)
    tᵣ::typeof(t) = solve(prob; reltol=1e-18, abstol=1e-18).u
    if tᵣ < zero(tᵣ)
        tᵣ = zero(tᵣ)
    end
    vtr = v(tᵣ, p)
    ztr = z(tᵣ, p)
    return ScalarField.Πs(x, ztr, vtr, p.q)
end

# function Π_ret_exact(x, p, t)
#     f(tr, p) = t - tr - abs(x - z(tr, p))
#     prob = NonlinearProblem(f, zero(t), p)
#     tr = solve(prob; reltol=1e-18, abstol=1e-18).u
#     tr = max(zero(tr), tr)

#     zr = z(tr, p)
#     vr = v(tr, p)
#     # g(xx, p) = -(p.q / 2) * sqrt(1 - vr^2) / (1 - sign(xx - zr) * vr)
#     # @show g(zr, p) - g(z(0.0, p), p)\
#     return -(p.q / 2) * sqrt(one(vr) - vr * vr) / (one(t) - sign(x - zr) * vr) #+ p.q / 2
# end

function Π_ret_exact(x, p, t)
    if t - abs(x - z(0.0, p)) <= 0
        return zero(t)
    end

    f(tr, p) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem(f, t, p)
    tr::typeof(t) = solve(prob; reltol=1e-18, abstol=1e-18).u

    zr = z(tr, p)
    vr = v(tr, p)
    s = sign(x - zr)

    return -0.5 * p.q * sqrt(one(vr) - vr * vr) /
           (one(vr) - s * vr)
end

function Ψ_ret(x, p, t)
    # calculate retarded time
    f(tr, p) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem(f, t, p)
    tᵣ::typeof(t) = solve(prob).u
    if tᵣ < zero(tᵣ)
        tᵣ = zero(tᵣ)
    end
    return ScalarField.Ψs(x, z(tᵣ, p), v(tᵣ, p), p.q)
end

# function Ψ_ret_exact(x, p, t)
#     # calculate retarded time
#     f(tr, p) = t - tr - abs(x - z(tr, p))
#     prob = NonlinearProblem(f, t, p)
#     tᵣ = solve(prob; reltol=1e-18, abstol=1e-18).u
#     if tᵣ < 0
#         tᵣ = zero(tᵣ)
#     end
#     @inline R(t, p) = (x - z(t, p))
#     return 0.5 * p.q *
#            (sign(R(tᵣ, p) * sqrt(one(v(tᵣ, p)) - v(tᵣ, p)^2)) /
#             (one(t) - v(tᵣ, p) * sign(R(tᵣ, p))))
# end

function Ψ_ret_exact(x, p, t)
    # First check whether the observation point is causally connected
    # to the source after switch-on at τ = 0.
    if t - abs(x - z(0.0, p)) <= 0
        return zero(t)
    end

    f(tr, p) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem(f, t, p)
    tr::typeof(t) = solve(prob; reltol=1e-18, abstol=1e-18).u

    zr = z(tr, p)
    vr = v(tr, p)
    s = sign(x - zr)

    return 0.5 * p.q * s * sqrt(one(vr) - vr * vr) /
           (one(vr) - s * vr)
end

end #end of module
