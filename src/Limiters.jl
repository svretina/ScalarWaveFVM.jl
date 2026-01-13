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

# @inline function superbee(σd, σu)
#     σ1 = minmod(σd, 2σu)
#     σ2 = minmod(2σd, σu)
#     return maxmod(σ1, σ2)
# end

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

@inline function unlimited(θ)
    return max(zero(θ), min((one(θ) + θ) * 0.5, 2.0, 2θ))
end

@inline function unlimited(a, b)
    return 0.5 * (a + b)
end

#### High resolution limiters:
@inline function minmod(θ)
    return minmod(one(θ), θ)
end

@inline function superbee(θ)
    return max(zero(θ), min(one(θ), 2θ), min(2.0, θ))
end

@inline function superbee(a, b)
    ab = a / b
    if isnan(ab)
        return zero(a)
    else
        return superbee(ab) * b
    end
end

@inline function MC(θ)
    return max(zero(θ), min((one(θ) + θ) * 0.5, 2.0, 2θ))
end

@inline function MC(a, b)
    ab = a / b
    if isnan(ab)
        return zero(a)
    else
        return MC(ab) * b
    end
end

# --- Van Leer ---
# Symmetric form: (r + |r|) / (1 + |r|)
# This correctly returns 0 for r <= 0, so the theta version was actually fine logically,
# but the (a,b) wrapper was unsafe.

@inline function vanLeer(θ)
    return (θ + abs(θ)) / (one(θ) + abs(θ))
end

@inline function vanLeer(a, b)
    # Optimized symmetric form: 2ab / (a+b) if signs match, else 0
    # checks if signs are same (a*b > 0) to avoid div by zero
    if a * b <= 0
        return zero(a)
    else
        return (2 * a * b) / (a + b)
    end
end

# --- Koren ---
# Koren is not symmetric, so we cannot easily simplify the (a,b) form without division.
# We must use the safety check.

@inline function koren(θ)
    return max(zero(θ), min(2θ, min((one(θ) + 2θ) / 3, 2.0)))
end

@inline function koren(a, b)
    # 1. Stability Check: If gradients oppose or b is 0, revert to 1st order (0 slope)
    if a * b <= 0
        return zero(a)
    end

    # 2. Koren Formula (The 3rd order sweet spot)
    # The term (a + 2b)/3 is what gives it the peak-preserving accuracy
    return min(2.0 * a, 2.0 * b, (a + 2.0 * b) / 3.0)
end

# --- Ospre ---
# FIXED: Added max(0, ...) to enforce TVD

@inline function ospre(θ)
    # Original: returns negative values for -1 < θ < 0
    # Fixed: Clamp to zero
    if θ <= 0
        return zero(θ)
    end
    θ2 = θ * θ
    return (1.5 * (θ2 + θ)) / (θ2 + θ + one(θ))
end

@inline function ospre(a, b)
    # Safety check
    if a * b <= 0
        return zero(a)
    end
    # We can pass through since we checked signs
    return ospre(a / b) * b
end

# --- Van Albada ---
# FIXED: Added max(0, ...) to enforce TVD

@inline function vanAlbada(θ)
    if θ <= 0
        return zero(θ)
    end
    θ2 = θ * θ
    return (θ2 + θ) / (θ2 + one(θ))
end

@inline function vanAlbada(a, b)
    # Symmetric form: (a^2*b + a*b^2) / (a^2 + b^2)
    # But usually simpler to just check signs
    if a * b <= 0
        return zero(a)
    end
    return vanAlbada(a / b) * b
end

@inline function umist(θ)
    # This formula is correct. It creates a "wide" second-order region.
    return max(zero(θ), min(2θ, 0.25 + 0.75θ, 0.75 + 0.25θ, 2.0))
end

@inline function umist(a, b)
    # FIX: Check for sign mismatch or b=0
    if a * b <= 0
        return zero(a)
    end
    return umist(a / b) * b
end

@inline function muscl(θ)
    return max(zero(θ), min(one(θ) + θ, 2.0, 2θ))
end

@inline function muscl(a, b)
    return muscl(a / b) * b
end

# "MP4-like" Limiter (Generalized Minmod with Theta=2)
# This mimics the compressive behavior of higher order schemes
@inline function generalized_minmod(a, b)
    theta = 2.0 # theta=1 is Minmod, theta=2 is Superbee/Roe-Compressive

    if a * b <= 0
        return zero(a)
    end

    # The "Mutidimensional" or "Generalized" limiter form
    # min(theta*|a|, |b|) ...

    # A very robust, sharp limiter often used as a substitute for full MP5
    # is the "Monotonized Central (MC)" you already have, OR "Superbee".

    # If you want something sharper than MC but safer than Superbee:
    return max(0.0, min(2.0 * abs(a), 2.0 * abs(b), 0.5 * (abs(a) + abs(b)))) * sign(a)
end

@inline function generalized_mc(a, b, theta::Float64=1.5)
    # 1. Stability Check
    if a * b <= 0
        return zero(a)
    end

    # 2. Generalized MC Formula
    # min( theta*|a|, theta*|b|, 0.5*|a+b| ) * sign(a)

    # We use the optimized form without expensive abs/sign calls where possible

    avg_slope = 0.5 * (a + b)

    # Apply the limiter bounds
    if b > 0
        return min(theta * a, theta * b, avg_slope)
    else
        return max(theta * a, theta * b, avg_slope)
    end
end

function create_tvb_limiter(base_limiter, M, h)
    # Pre-calculate the threshold to avoid doing M*h^2 every single cell
    # The threshold determines "How small must the change be to count as a smooth peak?"
    threshold = M * h^2

    return (a, b) -> begin
        # 1. TVB Logic:
        # If the difference between slopes is tiny, we are at a smooth extremum (peak/valley).
        # We bypass the limiter and return the Unlimited Central Difference.
        if abs(a - b) <= threshold
            return 0.5 * (a + b)
        end

        # 2. Fallback:
        # If the change is large (steep gradient or shock), use the strict TVD limiter.
        return base_limiter(a, b)
    end
end

end #end of module
