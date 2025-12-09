module RetardedPotentials

using NonlinearSolve
using ForwardDiff: derivative

export Φ_ret, Π_ret, Ψ_ret

@inline function retarded_time(tr, p)
    return -p.t + tr + abs(p.x - z(tr, p))
end

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
    tᵣ = solve(prob).u

    @inline R(t, p) = abs(x - z(t, p))
    # if tᵣ < 0
    #     return 0.5p.q * R(t, p)
    # else
    #     ϕ = 0.5p.q * (R(tᵣ, p) - v(tᵣ, p) * R(tᵣ, p))
    #     return ϕ
    # end

    ## retarded potential ?
    ϕ = 0.5p.q * (R(tᵣ, p) - v(tᵣ, p) * R(tᵣ, p)) /
        sqrt(1 - v(tᵣ, p)^2)
    ## singular field evaluated at retarded time
    ϕ = 0.5p.q * (R(tᵣ, p)) / sqrt(1 - v(tᵣ, p)^2)
    return ϕ
end

function Π_ret_autodiff(x, p, t)
    return derivative(tt -> Φ_ret(x, p, tt), t)
end

function Π_ret(x, p, t)
    # calculate retarded time
    f(tr, p) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem(f, 0.0, p)
    tᵣ = solve(prob).u

    @inline R(t, p) = (x - z(t, p))
    # if tᵣ < 0
    #     return -0.5 * p.q * v(t, p) * sign(R(t, p))
    # else
    #     return -0.5 * p.q * v(t, p) * sign(R(t, p)) -
    #            0.5 * p.q * v(tᵣ, p) * sign(R(tᵣ, p) - v(tᵣ, p) * R(tᵣ, p))
    # end

    ## time derivative of singular field evaluated at retarded time
    return -0.5 * p.q * v(tᵣ, p) * sign(R(tᵣ, p) - v(tᵣ, p) * R(tᵣ, p)) /
           sqrt(1 - v(tᵣ, p)^2)

    ## time derivative of retarded field
    #return 0.5 * p.q * (-v(tᵣ, p) + v(tᵣ, p)^2 - R(tᵣ, p) * a(tᵣ, p)) *
    #       sign(R(tᵣ, p) - v(tᵣ, p) * R(tᵣ, p))
    # return (v(tᵣ, p) - v(tᵣ, p)^2 + R(tᵣ, p) * a(tᵣ, p)) / ((1 - v(tᵣ, p))^3 * R(tᵣ, p)^2)
end

function Ψ_ret(x, p, t)
    # calculate retarded time
    f(tr, p) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem(f, 0.0, p)
    tᵣ = solve(prob).u

    @inline R(t, p) = (x - z(t, p))
    # if tᵣ < 0
    #     return -0.5 * p.q * v(t, p) * sign(R(t, p))
    # else
    #     return -0.5 * p.q * v(t, p) * sign(R(t, p)) -
    #            0.5 * p.q * v(tᵣ, p) * sign(R(tᵣ, p) - v(tᵣ, p) * R(tᵣ, p))
    # end
    return 0.5 * p.q * sign(R(tᵣ, p) - v(tᵣ, p) * R(tᵣ, p)) * (1 - v(tᵣ, p))
    # return 0.5p.q * sign(R(tᵣ, p) - v(tᵣ, p) * R(tᵣ, p))
end

end #end of module
