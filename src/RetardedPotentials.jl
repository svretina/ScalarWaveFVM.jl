module RetardedPotentials

using Roots: Bisection, find_zero
using QuadGK: quadgk

export Φ_ret, Π_ret, Ψ_ret, Φ_ret_exact, Π_ret_exact, Ψ_ret_exact

@inline z(t, p) = p.z0 + p.A * sin(p.ω * t)
@inline v(t, p) = p.ω * p.A * cos(p.ω * t)

"""
Return the coordinate retarded time for a source switched on at `t=0`.
`nothing` denotes an observation point outside the future light cone of the
initial worldline event.
"""
function retarded_time(x, p, t)
    t <= zero(t) && return nothing
    f(tr) = t - tr - abs(x - z(tr, p))
    f(zero(t)) <= zero(t) && return nothing
    return find_zero(f, (zero(t), t), Bisection())
end

function Φ_ret_exact(x, p, t)
    tr = retarded_time(x, p, t)
    tr === nothing && return zero(t)
    τr, _ = quadgk(s -> sqrt(one(s) - v(s, p)^2), zero(tr), tr)
    return -0.5 * p.q * τr
end

function Π_ret_exact(x, p, t)
    tr = retarded_time(x, p, t)
    tr === nothing && return zero(t)
    zr = z(tr, p)
    vr = v(tr, p)
    n = sign(x - zr)
    return -0.5 * p.q * sqrt(one(vr) - vr^2) / (one(vr) - n * vr)
end

function Ψ_ret_exact(x, p, t)
    tr = retarded_time(x, p, t)
    tr === nothing && return zero(t)
    zr = z(tr, p)
    vr = v(tr, p)
    n = sign(x - zr)
    return 0.5 * p.q * n * sqrt(one(vr) - vr^2) / (one(vr) - n * vr)
end

# Keep the historical public names, but make them exact aliases.  Evaluating
# the instantaneous puncture at retarded data is not the retarded potential.
Φ_ret(x, p, t) = Φ_ret_exact(x, p, t)
Π_ret(x, p, t) = Π_ret_exact(x, p, t)
Ψ_ret(x, p, t) = Ψ_ret_exact(x, p, t)

end
