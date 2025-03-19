module ScalarField

@inline function exp32(x::T) where {T<:Real}
    return x * sqrt(x)
end

# boost with +v

@inline function Φs(x, q, xp, v)
    return (q * abs(-x + xp)) / (2sqrt(1 - v^2))
end

@inline function Πs(x, q, xp, v)
    return -q * v * sign(x - xp) / (2sqrt(one(v) - v * v))
end

@inline function Ψs(x, q, xp, v)
    return q * sign(x - xp) / (2sqrt(one(v) - v * v))
end

@inline function ∂vΠ(x, q, xp, v)
    return -q * sign(x - xp) / (2exp32(one(v) - v * v))
end

@inline function ∂vΨ(x, q, xp, v)
    return v * q * sign(x - xp) / (2exp32(one(v) - v * v))
end

@inline function left_integral(q, xp, v, h)
    return h * q / (2exp32(one(v) - v * v))
end

@inline function left_cell_average(q, xp, v, h)
    return q / (2exp32(one(v) - v * v))
end

@inline function right_integral(q, xp, v, h)
    return -h * q / (2exp32(one(v) - v * v))
end

@inline function right_cell_average(q, xp, v, h)
    return -q / (2exp32(one(v) - v * v))
end

@inline function between_integral(left_face, q, xp, v, h)
    xi = left_face + 0.5h
    return q * (xp - xi) / (exp32(one(v) - v * v))
end

@inline function between_cell_average(left_face, q, xp, v, h)
    return between_integral(left_face, q, xp, v, h) / h
end

function cell_average(left_face, right_face, q, xp, v, h)
    if left_face < xp && right_face <= xp
        return left_cell_average(q, xp, v, h)
    elseif left_face >= xp && right_face > xp
        return right_cell_average(q, xp, v, h)
    elseif left_face < xp < right_face
        return between_cell_average(left_face, q, xp, v, h)
    else
        @show left_face
        @show right_face
        @show xp
        throw()
    end
    return nothing
end

# boost with -v
# @inline function Πs(x, q, xp, v)
#     return q * v * sign(x - xp) / (2sqrt(one(v) - v * v))
# end

# @inline function Ψs(x, q, xp, v)
#     return q * sign(x - xp) / (2sqrt(one(v) - v * v))
# end

# @inline function ∂vΠ(x, q, xp, v)
#     return q * sign(x - xp) / (2exp32(one(v) - v * v))
# end

# @inline function ∂vΨ(x, q, xp, v)
#     return v * q * sign(x - xp) / (2exp32(one(v) - v * v))
# end

end # end of module
