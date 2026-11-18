module ScalarField

using ForwardDiff: derivative

@inline function exp32(x::T) where {T<:Real}
    return x * sqrt(x)
end

# boost with +v

"""
    Φs(x, z, v, q)

Singular field convention used throughout the code, with `σ = sign(x - z)`:

`Φs = +(q/2) * |x - z| / (1 - vσ)`
`Πs = -(q/2) * vσ / (1 - vσ)`
`Ψs = +(q/2) * σ / (1 - vσ) = ∂ₓΦs`

The singular-source insertion in `ODE.add_singular_sources!` is derived from this
ansatz by adding the cell averages of `-a∂ᵥΠs`, `-a∂ᵥΨs`, and `-a∂ᵥΦs`.
"""
@inline function Φs(x, z, v, q, t)
    # return 0.5q * abs(x - z) * cosh(v)
    # return (q / 2) * abs(x - z) / (sqrt(one(v) - v * v))
    return 0.5 * q * sqrt(one(v) - v * v) * (abs(x - z) / (one(v) - v * sign(x - z)) - t)
end

@inline function Πs(x, z, v, q)
    # return -0.5q * sinh(v) * sign(x - z)
    # return -q * v * sign(x - z) / (2sqrt(one(v) - v * v))
    #  return -0.5 * q * v * sign(x - z) / (one(v) - v * sign(x - z))
    return -(q / 2) * sqrt(one(v) - v * v) / (one(v) - v * sign(x - z))
end

@inline function Ψs(x, z, v, q)
    return q / 2 * (sign(x - z) * sqrt(one(v) - v * v)) / (one(v) - sign(x - z) * v)
    # return 0.5 * q * sign(x - z) / (one(v) - v * sign(x - z))
    # return q * sign(x - z) / (2sqrt(one(v) - v * v))
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
    cell_average_Φs(xi, h, z, v, q)

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

@inline function cell_average_Φs(xi, h, z, v, q, t)
    left_face, right_face, _, _, lenL, lenR = _singular_cell_geometry(xi, h, z)

    if z >= right_face
        return 0.5 * (Φs(left_face, z, v, q, t) + Φs(right_face, z, v, q, t))
    elseif z <= left_face
        return 0.5 * (Φs(left_face, z, v, q, t) + Φs(right_face, z, v, q, t))
    else
        Φz = -0.5 * q * sqrt(one(v) - v * v) * t

        return 0.5 *
               (lenL * Φs(left_face, z, v, q, t) +
                lenR * Φs(right_face, z, v, q, t) +
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
@inline function cell_average_dvΦs(xi, h, z, v, q, t)
    return derivative(vv -> cell_average_Φs(xi, h, z, vv, q, t), v)
end

"""
    cell_average_dvΠs(xi, h, z, v, q)

Return the exact finite-volume cell average of `∂vΠs` over
`[xi - h/2, xi + h/2]`.
"""
@inline function cell_average_dvΠs(xi, h, z, v, q)
    return derivative(vv -> cell_average_Πs(xi, h, z, vv, q), v)
end

"""
    cell_average_dvΨs(xi, h, z, v, q)

Return the exact finite-volume cell average of `∂vΨs` over
`[xi - h/2, xi + h/2]`.
"""
@inline function cell_average_dvΨs(xi, h, z, v, q)
    return derivative(vv -> cell_average_Ψs(xi, h, z, vv, q), v)
end

"""
    cell_averages_Φs(x, h, z, v, q)

Return exact finite-volume cell averages of `Φs` over the cells centered at `x`.
"""
@inline cell_averages_Φs(x, h, z, v, q) = cell_average_Φs.(x, h, z, v, q)

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
    # return 0.5q * abs(x - z) * sinh(v) * cosh(v)^2
    return derivative(vv -> Φs(x, z, vv, q), v)
    # return 0.5q * v * abs(x - z) / exp32(one(v) - v * v)
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
