module RetardedPotentials

using NonlinearSolve
using Roots
using ForwardDiff: derivative
using ..ScalarField

export Φ_ret, Π_ret, Ψ_ret, Φ_ret2, Ψ_ret2

@inline function retarded_time(tr, p)
    return -p.t + tr + abs(p.x - z(tr, p))
end

@inline function z(t, p)
    # return 0.1t + 0.5 * 0.01 * t * t
    return p.z0 + p.A * sin(p.ω * t)
end

@inline function v(t, p)
    # return 0.1 + 0.01 * t
    return p.ω * p.A * cos(p.ω * t)
end

@inline function a(t, p)
    return -p.ω^2 * p.A * sin(p.ω * t)
end

function Φ_ret(x, p, t)
    # calculate retarded time
    f(tr, p) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem(f, zero(t), p)
    tᵣ = solve(prob; reltol=1e-12, abstol=1e-12).u
    if tᵣ < 0
        tᵣ = zero(tᵣ)
    end
    return ScalarField.Φs(x, z(tᵣ, p), v(tᵣ, p), p.q)
    # return -p.q * tᵣ / 2
end

function Π_ret(x, p, t)
    # calculate retarded time
    f(tr, p) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem(f, zero(t), p)
    tᵣ = solve(prob; reltol=1e-18, abstol=1e-18).u
    # @show tᵣ
    if tᵣ < 0
        tᵣ = zero(tᵣ)
    end
    @inline R(t, p) = x - z(t, p)
    return ScalarField.Πs(x, z(tᵣ, p), v(tᵣ, p), p.q)
    # return -0.5p.q * (1 / (1 - v(tᵣ, p) * sign(R(tᵣ, p)))) + 0.5
end

function Ψ_ret(x, p, t)
    # calculate retarded time
    f(tr, p) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem(f, t, p)
    tᵣ = solve(prob).u
    if tᵣ < 0
        tᵣ = zero(tᵣ)
    end
    @inline R(t, p) = (x - z(t, p))
    return ScalarField.Ψs(x, z(tᵣ, p), v(tᵣ, p), p.q)
    # return 0.5p.q * (sign(R(tᵣ, p)) / (one(t) - v(tᵣ, p) * sign(R(tᵣ, p))))
end

end #end of module
