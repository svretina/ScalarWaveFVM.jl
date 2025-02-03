module Limiters

using ..Equations

using StaticArrays

@inline function minmod(σ1::T, σ2::T) where {T}
    if sign(σ1) == sign(σ2)
        return sign(σ1) * min(abs(σ1), abs(σ2))
    else
        return zero(σ1)
    end
end

# @inline function minmod(σ1::T, σ2::T) where {T}
#     if σ1 > 0 && σ2 > 0
#         return min(σ1, σ2)
#     elseif σ1 < 0 && σ2 < 0
#         return max(σ1, σ2)
#     else
#         return zero(σ1)
#     end
# end

function minmod(a::T, b::T, c::T) where {T}
    if sign(a) == sign(b) == sign(c)
        return sign(a) * min(abs(a), abs(b), abs(c))
    else
        return zero(a)
    end
end

@inline function maxmod(σ1, σ2)
    if sign(σ1) == sign(σ2)
        return sign(σ1) * max(abs(σ1), abs(σ2))
    else
        return zero(σ1)
    end
end

## Slope limiters

@inline function superbee(σd, σu)
    σ1 = minmod(σd, 2σu)
    σ2 = minmod(2σd, σu)
    return maxmod(σ1, σ2)
end

@inline function upwind(σ1::T, σ2::T) where {T}
    return zero(σ1)
end

### Flux Limiters

### Linear methods limiters:
@inline function upwind(θ)
    return zero(θ)
end

@inline function Lax_Wendroff(θ)
    return one(θ)
end

@inline function Beam_Warming(θ)
    return θ
end

@inline function Fromm(θ)
    return 0.5 * (one(θ) + θ)
end

#### High resolution limiters:
@inline function minmod(θ)
    return minmod(one(θ), θ)
end

@inline function superbee(θ)
    return max(zero(θ), min(one(θ), 2θ), min(2.0, θ))
end

@inline function MC(θ)
    return max(zero(θ), min((one(θ) + θ) * 0.5, 2.0, 2θ))
end

@inline function vanLeer(θ)
    return (θ + abs(θ)) / (one(θ) + abs(θ))
end

@inline function koren(θ)
    return max(zero(θ), min(2θ, min((one(θ) + 2θ) / 3, 2.0)))
end

@inline function ospre(θ)
    θ2 = θ * θ
    return (1.5 * (θ2 + θ)) / (θ2 + θ + one(θ))
end

@inline function vanAlbada(θ)
    θ2 = θ * θ
    return (θ2 + θ) / (θ2 + one(θ))
end

@inline function umist(θ)
    return max(zero(θ), min(2θ, 0.25 + 0.75θ, 0.75 + 0.25θ, 2.0))
end

@inline function muscl(θ)
    return max(zero(θ), min(one(θ) + θ, 2.0, 2θ))
end

end #end of module
