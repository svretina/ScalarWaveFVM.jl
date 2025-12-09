module ScalarField

using ForwardDiff

@inline function exp32(x::T) where {T<:Real}
    return x * sqrt(x)
end

# boost with +v

@inline function Φs(x, z, v, q)
    #return 0.5q * abs(x - z) * cosh(v)
    return (q / 2) * abs((x - z) / (sqrt(one(v) - v * v)))
end

@inline function Πs(x, z, v, q)
    # return -0.5q * sinh(v) * sign(x - z)
    return -q * v * sign(x - z) / (2sqrt(one(v) - v * v))
end

@inline function Ψs(x, z, v, q)
    return q * sign(x - z) / (2sqrt(one(v) - v * v))
end

## Φ
@inline function Φs001(x, z, v, q)
    # return 0.5q * abs(x - z) * sinh(v) * cosh(v)^2
    return derivative(vv -> Φs(x, z, vv, q), v)
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

## Π = ∂ₜΦ
@inline function Πs001(x, z, v, q)
    # return -0.5q * cosh(v)^3 * sign(x - z)
    return derivative(vv -> Πs(x, z, vv, q), v)
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

@inline function ∂²vΠs2(x, z, v, q)
    return -3q * v * sign(x - z) / (2(one(v) - v * v)^2.5)
end

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
