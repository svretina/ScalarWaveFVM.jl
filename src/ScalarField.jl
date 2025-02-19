module ScalarField

@inline function exp32(x::T) where {T<:Real}
    return x * sqrt(x)
end

function ∂vΠ(x, q, xp, v)
    return -q * sign(x - xp) / (2exp32(one(v) - v * v))
end

function ∂vΨ(x, q, xp, v)
    return v * q * sign(x - xp) / (2exp32(one(v) - v * v))
end

end # end of module
