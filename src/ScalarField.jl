module ScalarField

using ForwardDiff: derivative

@inline function exp32(x::T) where {T<:Real}
    return x * sqrt(x)
end

"""
    Φs(x, z, v, q, τ)

Retarded-matched 1+1 puncture for the positive-source convention
`□Φ = q δ(x-z)/γ`, with `r = x-z` and accumulated proper time `τ`:

`Φs = (qγ/2) * (|r| + v*r) - (q/2) * τ`
`Πs = -(qγ/2) * (1 + v*sign(r))`
`Ψs = +(qγ/2) * (sign(r) + v)`

The first-order fields equal the one-sided local retarded fields. Velocity
derivatives hold `r` and `τ` fixed, as required by the remainder equations.
"""
@inline function Φs(x, z, v, q, τ)
    r = x - z
    γ = inv(sqrt(one(v) - v * v))
    return 0.5 * q * (γ * (abs(r) + v * r) - τ)
end

@inline function Πs(x, z, v, q)
    n = sign(x - z)
    γ = inv(sqrt(one(v) - v * v))
    return -0.5 * q * γ * (one(v) + v * n)
end

@inline function Ψs(x, z, v, q)
    n = sign(x - z)
    γ = inv(sqrt(one(v) - v * v))
    return 0.5 * q * γ * (n + v)
end

@inline function dvΦs(x, z, v, q)
    r = x - z
    γ = inv(sqrt(one(v) - v * v))
    return 0.5 * q * γ^3 * (r + v * abs(r))
end

@inline function dvΠs(x, z, v, q)
    n = sign(x - z)
    γ = inv(sqrt(one(v) - v * v))
    return -0.5 * q * γ^3 * (v + n)
end


@inline function dvΨs(x, z, v, q)
    n = sign(x - z)
    γ = inv(sqrt(one(v) - v * v))
    return 0.5 * q * γ^3 * (one(v) + v * n)
end

"""Return `(Pi_s, Psi_s)` at a face, selecting `side` only when `x == z`."""
@inline function singular_face_trace(x, z, v, q, side)
    n = x < z ? -one(v) : x > z ? one(v) : sign(side)
    γ = inv(sqrt(one(v) - v * v))
    return -0.5 * q * γ * (one(v) + v * n),
           0.5 * q * γ * (n + v)
end

"""Return the left/right one-sided traces `(PiL, PiR, PsiL, PsiR)`."""
@inline function singular_particle_traces(z, v, q)
    ΠL, ΨL = singular_face_trace(z, z, v, q, -one(v))
    ΠR, ΨR = singular_face_trace(z, z, v, q, one(v))
    return ΠL, ΠR, ΨL, ΨR
end

"""
    initialize_remainder(x; mode, z, v, q, τ, full_field, custom)

Construct cell-centred remainder data in `(Πr, Ψr, Φr)` column order.
`mode=:zero_remainder` selects the global puncture extension, `:retarded`
subtracts the puncture from a supplied full-field callback, and `:custom`
copies user-supplied data. The callback must return `(Φ, Π, Ψ)`.
"""
function initialize_remainder(x; mode=:zero_remainder, z=0.0, v=0.0, q=0.0,
                              τ=0.0, full_field=nothing, custom=nothing)
    if mode === :zero_remainder
        return zeros(promote_type(eltype(x), typeof(q)), length(x), 3)
    elseif mode === :custom
        custom === nothing && throw(ArgumentError("custom remainder data are required"))
        size(custom) == (length(x), 3) ||
            throw(DimensionMismatch("custom data must have size ($(length(x)), 3)"))
        return copy(custom)
    elseif mode === :retarded
        full_field === nothing &&
            throw(ArgumentError("full_field callback is required for :retarded initialization"))
        data = Matrix{promote_type(eltype(x), typeof(q))}(undef, length(x), 3)
        for (i, xi) in pairs(x)
            full = full_field(xi)
            data[i, 1] = full[2] - Πs(xi, z, v, q)
            data[i, 2] = full[3] - Ψs(xi, z, v, q)
            data[i, 3] = full[1] - Φs(xi, z, v, q, τ)
        end
        return data
    end
    throw(ArgumentError("unknown initialization mode $mode"))
end

@inline function _cell_faces(xi, h)
    half_h = 0.5 * h
    return xi - half_h, xi + half_h
end

@inline function _singular_cell_geometry(xi, h, z)
    left_face, right_face = _cell_faces(xi, h)
    dxL = z - left_face
    dxR = right_face - z
    lenL = clamp(dxL, zero(h), h)
    lenR = h - lenL
    return left_face, right_face, dxL, dxR, lenL, lenR
end

"""
    cell_average_Φs(xi, h, z, v, q, τ)

Return the exact finite-volume cell average of `Φs` over
`[xi - h/2, xi + h/2]`.
"""
# @inline function cell_average_Φs(xi, h, z, v, q, t)
#     left_face, right_face, _, _, lenL, lenR = _singular_cell_geometry(xi, h, z)

#     if z >= right_face
#         return 0.5 * (Φs(left_face, z, v, q, t) + Φs(right_face, z, v, q, t))
#     elseif z <= left_face
#         return 0.5 * (Φs(left_face, z, v, q, t) + Φs(right_face, z, v, q, t))
#     else
#         return 0.5 *
#                (lenL * Φs(left_face, z, v, q, t) + lenR * Φs(right_face, z, v, q, t)) / h
#     end
# end

@inline function cell_average_Φs(xi, h, z, v, q, τ)
    left_face, right_face, _, _, lenL, lenR = _singular_cell_geometry(xi, h, z)

    if z >= right_face
        return 0.5 * (Φs(left_face, z, v, q, τ) + Φs(right_face, z, v, q, τ))
    elseif z <= left_face
        return 0.5 * (Φs(left_face, z, v, q, τ) + Φs(right_face, z, v, q, τ))
    else
        Φz = -0.5 * q * τ

        return 0.5 *
               (lenL * Φs(left_face, z, v, q, τ) +
                lenR * Φs(right_face, z, v, q, τ) +
                h * Φz) / h
    end
end

"""
    cell_average_Πs(xi, h, z, v, q)

Return the exact finite-volume cell average of `Πs` over
`[xi - h/2, xi + h/2]`.
"""
@inline function cell_average_Πs(xi, h, z, v, q)
    left_face, right_face, _, _, lenL, lenR = _singular_cell_geometry(xi, h, z)

    if z >= right_face
        return Πs(xi, z, v, q)
    elseif z <= left_face
        return Πs(xi, z, v, q)
    else
        xL = 0.5 * (left_face + z)
        xR = 0.5 * (z + right_face)
        return (lenL * Πs(xL, z, v, q) + lenR * Πs(xR, z, v, q)) / h
    end
end

"""
    cell_average_Ψs(xi, h, z, v, q)

Return the exact finite-volume cell average of `Ψs` over
`[xi - h/2, xi + h/2]`.
"""
@inline function cell_average_Ψs(xi, h, z, v, q)
    left_face, right_face, _, _, lenL, lenR = _singular_cell_geometry(xi, h, z)

    if z >= right_face
        return Ψs(xi, z, v, q)
    elseif z <= left_face
        return Ψs(xi, z, v, q)
    else
        xL = 0.5 * (left_face + z)
        xR = 0.5 * (z + right_face)
        return (lenL * Ψs(xL, z, v, q) + lenR * Ψs(xR, z, v, q)) / h
    end
end

"""
    cell_average_dvΦs(xi, h, z, v, q)

Return the exact finite-volume cell average of `∂vΦs` over
`[xi - h/2, xi + h/2]`.
"""
@inline function cell_average_dvΦs(xi, h, z, v, q)
    left_face, right_face, _, _, lenL, lenR = _singular_cell_geometry(xi, h, z)
    if z >= right_face || z <= left_face
        return dvΦs(xi, z, v, q)
    end
    return 0.5 * (lenL * dvΦs(left_face, z, v, q) +
                  lenR * dvΦs(right_face, z, v, q)) / h
end

"""
    cell_average_dvΠs(xi, h, z, v, q)

Return the exact finite-volume cell average of `∂vΠs` over
`[xi - h/2, xi + h/2]`.
"""
@inline function cell_average_dvΠs(xi, h, z, v, q)
    _, _, _, _, lenL, lenR = _singular_cell_geometry(xi, h, z)
    γ = inv(sqrt(one(v) - v * v))
    left = -0.5 * q * γ^3 * (v - one(v))
    right = -0.5 * q * γ^3 * (v + one(v))
    return (lenL * left + lenR * right) / h
end

"""
    cell_average_dvΨs(xi, h, z, v, q)

Return the exact finite-volume cell average of `∂vΨs` over
`[xi - h/2, xi + h/2]`.
"""
@inline function cell_average_dvΨs(xi, h, z, v, q)
    _, _, _, _, lenL, lenR = _singular_cell_geometry(xi, h, z)
    γ = inv(sqrt(one(v) - v * v))
    left = 0.5 * q * γ^3 * (one(v) - v)
    right = 0.5 * q * γ^3 * (one(v) + v)
    return (lenL * left + lenR * right) / h
end

"""
    cell_averages_Φs(x, h, z, v, q, τ)

Return exact finite-volume cell averages of `Φs` over the cells centered at `x`.
"""
@inline cell_averages_Φs(x, h, z, v, q, τ) = cell_average_Φs.(x, h, z, v, q, τ)

"""
    cell_averages_Πs(x, h, z, v, q)

Return exact finite-volume cell averages of `Πs` over the cells centered at `x`.
"""
@inline cell_averages_Πs(x, h, z, v, q) = cell_average_Πs.(x, h, z, v, q)

"""
    cell_averages_Ψs(x, h, z, v, q)

Return exact finite-volume cell averages of `Ψs` over the cells centered at `x`.
"""
@inline cell_averages_Ψs(x, h, z, v, q) = cell_average_Ψs.(x, h, z, v, q)

"""
    cell_averages_dvΦs(x, h, z, v, q)

Return exact finite-volume cell averages of `∂vΦs` over the cells centered at `x`.
"""
@inline cell_averages_dvΦs(x, h, z, v, q) = cell_average_dvΦs.(x, h, z, v, q)

"""
    cell_averages_dvΠs(x, h, z, v, q)

Return exact finite-volume cell averages of `∂vΠs` over the cells centered at `x`.
"""
@inline cell_averages_dvΠs(x, h, z, v, q) = cell_average_dvΠs.(x, h, z, v, q)

"""
    cell_averages_dvΨs(x, h, z, v, q)

Return exact finite-volume cell averages of `∂vΨs` over the cells centered at `x`.
"""
@inline cell_averages_dvΨs(x, h, z, v, q) = cell_average_dvΨs.(x, h, z, v, q)

## Φ
@inline function Φs001(x, z, v, q)
    return dvΦs(x, z, v, q)
end

@inline function Φs002(x, z, v, q)
    return derivative(vv -> Φs001(x, z, vv, q), v)
end

## Ψ = ∂ₓΦ
@inline function Ψs100(x, z, v, q)
    return derivative(xx -> Ψs(xx, z, v, q), x)
end
@inline function Ψs101(x, z, v, q)
    return derivative(vv -> Ψs100(x, z, vv, q), v)
end
@inline function Ψs001(x, z, v, q)
    return derivative(vv -> Ψs(x, z, vv, q), v)
end

@inline function Ψs002(x, z, v, q)
    return derivative(vv -> Ψs001(x, z, vv, q), v)
end

@inline function Πs001(x, z, v, q)
    # return -0.5q * cosh(v)^3 * sign(x - z)
    return derivative(vv -> Πs(x, z, vv, q), v)
    # return -0.5q * sign(x - z) / exp32(one(v) - v * v)
end

@inline function Πs010(x, z, v, q)
    return derivative(zz -> Πs(x, zz, v, q), z)
end

@inline function Πs011(x, z, v, q)
    return derivative(vv -> Πs010(x, z, vv, q), v)
end

@inline function Πs002(x, z, v, q)
    return derivative(vv -> Πs001(x, z, vv, q), v)
end

# @inline function ∂xΨs(x, z, v, q)
#     return derivative(xx -> Ψs(xx, z, v, q), x)
# end

# @inline function ∂vΦs(x, z, v, q)
#     return derivative(vv -> Φs(x, z, vv, q), v)
# end

# @inline function ∂vΠs(x, z, v, q)
#     return derivative(vv -> Πs(x, z, vv, q), v)
# end

# @inline function ∂zΠs(x, z, v, q)
#     return derivative(zz -> Πs(x, zz, v, q), z)
# end

# @inline function ∂²vΠs(x, z, v, q)
#     return derivative(vv -> ∂vΠs(x, z, vv, q), v)
# end

# @inline function ∂vΨs(x, z, v, q)
#     return derivative(vv -> Ψs(x, z, vv, q), v)
# end

#####
# @inline function ∂xΨs2(x, z, v, q)
#     delta(x, z) = x == z ? one(x) : zero(x)
#     return q * delta(x, z) / (sqrt(one(v) - v * v))
# end

# @inline function ∂vΦs2(x, z, v, q)
#     return q * v * abs(x - z) / (2exp32(one(v) - v * v))
# end

# @inline function ∂vΠs2(x, z, v, q)
#     return -q * sign(x - z) / (2exp32(one(v) - v * v))
# end

# @inline function ∂zΠs2(x, z, v, q)
#     delta(x, z) = x == z ? one(x) : zero(x)
#     return q * v * delta(x, z) / (sqrt(one(v) - v * v))
# end

# @inline function ∂²vΠs2(x, z, v, q)
#     return -3q * v * sign(x - z) / (2(one(v) - v * v)^2.5)
# end

# @inline function ∂vΨs2(x, z, v, q)
#     return v * q * sign(x - z) / (2exp32(one(v) - v * v))
# end

# ## integrals and averages of ∂vΠ
# @inline function left_integral(q, z, v, h)
#     return h * q / (2exp32(one(v) - v * v))
# end

# @inline function left_cell_average(q, z, v, h)
#     return q / (2exp32(one(v) - v * v))
# end

# @inline function right_integral(q, z, v, h)
#     return -h * q / (2exp32(one(v) - v * v))
# end

# @inline function right_cell_average(q, z, v, h)
#     return -q / (2exp32(one(v) - v * v))
# end

# @inline function between_integral(left_face, q, z, v, h)
#     xi = left_face + 0.5h
#     return q * (z - xi) / (exp32(one(v) - v * v))
# end

# @inline function between_cell_average(left_face, q, z, v, h)
#     return between_integral(left_face, q, z, v, h) / h
# end

# function cell_average(left_face, right_face, q, z, v, h)
#     if left_face < z && right_face <= z
#         return left_cell_average(q, z, v, h)
#     elseif left_face >= z && right_face > z
#         return right_cell_average(q, z, v, h)
#     elseif left_face < z < right_face
#         return between_cell_average(left_face, q, z, v, h)
#     else
#         @show left_face
#         @show right_face
#         @show z
#         throw()
#     end
#     return nothing
# end

# boost with -v
# @inline function Πs(x, z, v, q)
#     return q * v * sign(x - z) / (2sqrt(one(v) - v * v))
# end

# @inline function Ψs(x, z, v, q)
#     return q * sign(x - z) / (2sqrt(one(v) - v * v))
# end

# @inline function ∂vΠ(x, z, v, q)
#     return q * sign(x - z) / (2exp32(one(v) - v * v))
# end

# @inline function ∂vΨ(x, z, v, q)
#     return v * q * sign(x - z) / (2exp32(one(v) - v * v))
# end

end # end of module
