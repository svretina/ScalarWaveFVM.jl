module InitialData

@inline function SineWave(t, x, A, λ, c)
    return A * sin(2π * (x - c * t) / λ)
end

@inline function dtSineWave(t::Real, x::Real, A, λ, c)
    return -2A * π * c * cos(2π * (x - c * t) / λ) / λ
end

@inline function dxSineWave(t::Real, x::Real, A, λ, c)
    return 2A * π * cos(2π * (x - c * t) / λ) / λ
end

@inline function Harmonic(x, k)
    return 0.5k * x * x
end

@inline function dxHarmonic(x, k)
    return k * x
end

@inline function sign(x, z)
    if x < z
        return one(x)
    elseif x > z
        return -one(x)
    else
        return zero(x)
    end
end

# @inline function Gaussian1D(t::Real, x::Real, A::Real, σ::Real)
#     L = 100
#     xc = t
#     xc = mod(xc + L, 2L) - L

#     x_shifted1 = x - xc + 2L
#     g1 = A * exp(-x_shifted1^2 / (2 * σ^2))
#     x_shiftedm1 = x - xc - 2L
#     gm1 = A * exp(-x_shiftedm1^2 / (2 * σ^2))
#     x_shifted0 = x - xc
#     g0 = A * exp(-x_shifted0^2 / (2 * σ^2))
#     return g1 + gm1 + g0
# end

# @inline function dtGaussian1D(t::Real, x::Real, A::Real, σ::Real)
#     return ForwardDiff.derivative(t1 -> Gaussian1D(t1, x, A, σ), t)
# end

# @inline function dxGaussian1D(t::Real, x::Real, A::Real, σ::Real)
#     return ForwardDiff.derivative(x1 -> Gaussian1D(t, x1, A, σ), x)
# end

end # end of module
