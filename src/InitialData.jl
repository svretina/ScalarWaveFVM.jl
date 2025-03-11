module InitialData

using ForwardDiff

function triangular_pulse(x, t, L, amplitude=1.0, base_width=2.0, speed=1.0)
    # x: spatial coordinate
    # t: time
    # L: half-length of the domain (domain is [-L, L])
    # amplitude: height of the triangle
    # base_width: width of the triangle base
    # speed: velocity of the pulse

    # Calculate the domain length
    domain_length = 2L

    # Calculate the position of the triangle center
    center = mod(speed * t - L, domain_length) - L

    # Calculate the left and right edges of the triangle
    left_edge = center - base_width / 2
    right_edge = center + base_width / 2

    # Adjust x to be within the domain
    x_adjusted = mod(x + L, domain_length) - L

    # Check if x is within the triangle, considering periodic boundaries
    if left_edge <= x_adjusted <= right_edge
        # Calculate the height at this point
        if x_adjusted <= center
            return amplitude * (x_adjusted - left_edge) / (base_width / 2)
        else
            return amplitude * (right_edge - x_adjusted) / (base_width / 2)
        end
    elseif left_edge < -L && x_adjusted >= left_edge + domain_length
        # Handle wrapping on the left side
        return amplitude * (x_adjusted - (left_edge + domain_length)) / (base_width / 2)
    elseif right_edge > L && x_adjusted <= right_edge - domain_length
        # Handle wrapping on the right side
        return amplitude * ((right_edge - domain_length) - x_adjusted) / (base_width / 2)
    else
        return 0.0
    end
end

function dt_triangular_pulse(x, t, L, amplitude=1.0, pulse_width=1.0, speed=1.0)
    return ForwardDiff.derivative(t1 -> triangular_pulse(x, t1, L, amplitude,
                                                         pulse_width, speed), t)
end

function dx_triangular_pulse(x, t, L, amplitude=1.0, pulse_width=1.0, speed=1.0)
    return ForwardDiff.derivative(x1 -> triangular_pulse(x1, t, L, amplitude,
                                                         pulse_width, speed), x)
end

function square_pulse(x, t, L, amplitude=1.0, pulse_width=1.0, speed=1.0)
    # x: spatial coordinate
    # t: time
    # L: half-length of the domain (domain is [-L, L])
    # amplitude: height of the pulse
    # pulse_width: width of the pulse
    # speed: velocity of the pulse

    # Calculate the domain length
    domain_length = 2L

    # Calculate the position of the pulse center
    center = mod(speed * t - L, domain_length) - L

    # Calculate the left and right edges of the pulse
    left_edge = center - pulse_width / 2
    right_edge = center + pulse_width / 2

    # Adjust x to be within the domain
    x_adjusted = mod(x + L, domain_length) - L

    # Check if x is within the pulse, considering periodic boundaries
    if (left_edge <= x_adjusted <= right_edge) ||
       (left_edge < -L && x_adjusted >= left_edge + domain_length) ||
       (right_edge > L && x_adjusted <= right_edge - domain_length)
        return amplitude
    else
        return 0.0
    end
end

function dtSquarePulse(x, t, L, amplitude=1.0, pulse_width=1.0, speed=1.0)
    return ForwardDiff.derivative(t1 -> square_pulse(x, t1, L, amplitude,
                                                     pulse_width, speed), t)
end

function dxSquarePulse(x, t, L, amplitude=1.0, pulse_width=1.0, speed=1.0)
    return ForwardDiff.derivative(x1 -> square_pulse(x1, t, L, amplitude,
                                                     pulse_width, speed), x)
end

@inline function Gaussian1D(t::Real, x::Real, A::Real, σ::Real)
    L = 100
    xc = t
    xc = mod(xc + L, 2L) - L

    x_shifted1 = x - xc + 2L
    g1 = A * exp(-x_shifted1^2 / (2 * σ^2))
    x_shiftedm1 = x - xc - 2L
    gm1 = A * exp(-x_shiftedm1^2 / (2 * σ^2))
    x_shifted0 = x - xc
    g0 = A * exp(-x_shifted0^2 / (2 * σ^2))
    return g1 + gm1 + g0
end

@inline function dtGaussian1D(t::Real, x::Real, A::Real, σ::Real)
    return ForwardDiff.derivative(t1 -> Gaussian1D(t1, x, A, σ), t)
end

@inline function dxGaussian1D(t::Real, x::Real, A::Real, σ::Real)
    return ForwardDiff.derivative(x1 -> Gaussian1D(t, x1, A, σ), x)
end

@inline function SineWave(t, x, A, λ, c)
    return A * sin(2π * (x - c * t) / λ)
end

@inline function dtSineWave(t::Real, x::Real, A, λ, c)
    return -2A * π * c * cos(2π * (x - c * t) / λ) / λ
end

@inline function dxSineWave(t::Real, x::Real, A, λ, c)
    return 2A * π * cos(2π * (x - c * t) / λ) / λ
end

@inline function Harmonic(t, x)
    return 0.5x * x
end

@inline function dtHarmonic(t, x)
    return zero(x)
end

@inline function dxHarmonic(t, x, a)
    return a * x
end

end # end of module
