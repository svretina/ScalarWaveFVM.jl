module Limiters

using ..Equations

using StaticArrays

## Linear methods limiters:
##
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

# measure of data smoothness
# θ ≈ 1 if smooth
# else far from 1
function θ(i, Q, equation::LinearScalarWaveEquation1D, direction::Int)
    if direction > 0
        I = i - 1
    elseif direction < 0
        I = i + 1
    else
        throw("direction cannot be $direction . Can be either -1 or +1")
    end
    Δqf = @SVector [Q[I, 1] - Q[I - 1, 1], Q[I, 2] - Q[I - 1, 2]]
    Δq = @SVector [Q[i, 1] - Q[i - 1, 1], Q[i, 2] - Q[i - 1, 2]]

    if Δq[1] == Δqf[1] && Δq[2] == Δqf[2]
        return @SVector [1.0, 1.0]
    elseif Δq[1] == Δqf[1] && Δq[2] !== Δqf[2]
        return @SVector [1.0, Δqf[2] / Δq[2]]
    elseif Δq[1] !== Δqf[1] && Δq[2] == Δqf[2]
        return @SVector [Δqf[1] / Δq[1], 1.0]
    else
        return Δqf ./ Δq
    end
end

## Slope limiters
function minmod(u, i, h)
    upwind = forward_slope(u, i, h)
    downwind = backward_slope(u, i, h)
    return minmod(upwind, downwind)
end

function minmod(σ1, σ2)
    if sign(σ1) == sign(σ2)
        return min(σ1, σ2)
    else
        return zero(σ1)
    end
end

function maxmod(σ1, σ2)
    if sign(σ1) == sign(σ2)
        return max(σ1, σ2m)
    else
        return zero(σ1)
    end
end

function superbee(u, i, h)
    σ1 = backward_slope(u, i, h)
    σ2 = forward_slope(u, i, h)
    σL = minmod(2σ1, σ2)
    σR = minmod(σ1, 2σ2)
    return maxmod(σL, σR)
end

## Flux Limiters

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
    return (θ + abs(θ)) / (one(θ) + θ)
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
