module Energies

using ..ScalarField
using ..ParticleMotion
using DataInterpolations

function newton_cotes(y::Vector{T}, h::T; order::Int=1) where {T}
    n = length(y) - 1  # Number of intervals
    if n < order
        error("Number of points ($(n+1)) must be at least $(order+1) for order $order Newton-Cotes.")
    end
    if (n % order) != 0
        error("Number of intervals ($n) must be a multiple of $order for order $order Newton-Cotes.")
    end

    integral = zero(T)

    if order == 1  # Trapezoidal Rule
        @inbounds for i in 1:n
            integral += y[i]
        end
        @inbounds integral += 0.5 * (y[1] + y[n + 1])  # Adjust endpoints
        integral *= h

    elseif order == 2  # Simpson's Rule
        @inbounds for i in 1:(n ÷ 2)
            j = 2 * i - 1  # Index of odd terms (4 coefficient)
            integral += 4.0 * y[j + 1]  # Middle point
            integral += y[j]            # Even point (coefficient 1 or 2 from next iteration)
        end
        @inbounds integral += y[n + 1]  # Add the last point (coefficient 1)
        integral *= h / 3.0

    elseif order == 3  # Simpson's 3/8 Rule
        @inbounds for i in 0:(n ÷ 3 - 1)
            j = 3 * i + 1
            integral += y[j]          # Coefficient 1
            integral += 3.0 * y[j + 1]  # Coefficient 3
            integral += 3.0 * y[j + 2]  # Coefficient 3
            integral += y[j + 3]      # Coefficient 1
        end
        integral *= 3.0 * h / 8.0

    else
        error("Order $order not implemented. Supported orders: 1, 2, 3.")
    end
    return integral
end

function HamiltonianDensity(i, sim)
    if hasproperty(sim.sol.u[1], :x)
        Ufield = sim.sol.u[i].x[1]
        Upart = sim.sol.u[i].x[2]
        x = sim.params.x
        q1 = sim.params.q1
        q2 = sim.params.q2

        xp1 = Upart[2]
        xp2 = Upart[5]
        vp1 = Upart[3]
        vp2 = Upart[6]

        Φs1 = ScalarField.Φs.(x, q1, xp1, vp1)
        Πs1 = ScalarField.Πs.(x, q1, xp1, vp1)
        Ψs1 = ScalarField.Ψs.(x, q1, xp1, vp1)

        Φs2 = ScalarField.Φs.(x, q2, xp2, vp2)
        Πs2 = ScalarField.Πs.(x, q2, xp2, vp2)
        Ψs2 = ScalarField.Ψs.(x, q2, xp2, vp2)
        Πr = @views Ufield[:, 1]
        Ψr = @views Ufield[:, 2]
        p = @. -(Πr + Πs1 + Πs2)
        H = @. -0.5p * p - p * (Πs1 + Πs2) - 0.5Ψr * Ψr - Ψr * Ψs1 - Ψr * Ψs2 -
               0.5Ψs1 * Ψs1 - 0.5Ψs2 * Ψs2 - Ψs1 * Ψs2
        return H
    else
        Ufield = sim.sol.u[i]
        x = sim.params.x
        q = sim.params.q
        xp, v, _ = ParticleMotion.oscillator(sim.sol.t[i], sim.params.x0, sim.params.A,
                                             sim.params.ω)
        # Φs = ScalarField.Φs.(x, q, xp, v)
        Πs = ScalarField.Πs.(x, q, xp, v)
        Ψs = ScalarField.Ψs.(x, q, xp, v)

        Πr = @views Ufield[:, 1]
        Ψr = @views Ufield[:, 2]

        p = -(Πr .+ Πs)
        H = @. -0.5p * p - p * Πs - 0.5Ψr * Ψr - Ψr * Ψs - 0.5Ψs * Ψs
        return H
    end
end

function Hamiltonian(i, sim)
    H = HamiltonianDensity(i, sim)
    energy = newton_cotes(H, sim.params.dx; order=1)

    if hasproperty(sim.sol.u[1], :x)
        Ufield = sim.sol.u[i].x[1]
        Upart = sim.sol.u[i].x[2]
        x = sim.params.x
        Φr = @views Ufield[:, 3]

        q1 = sim.params.q1
        q2 = sim.params.q2

        xp1 = Upart[2]
        xp2 = Upart[5]
        vp1 = Upart[3]
        vp2 = Upart[6]

        Φs1_at_2 = ScalarField.Φs.(xp2, q1, xp1, vp1)
        Φs2_at_1 = ScalarField.Φs.(xp1, q2, xp2, vp2)

        interpolator_Φr = QuadraticInterpolation(Φr, x)
        Φrz1 = interpolator_Φr(xp1)
        Φrz2 = interpolator_Φr(xp2)

        energy += -q1 * Φrz1 - q2 * Φrz2 - q1 * Φs2_at_1 - q2 * Φs1_at_2
        return energy
    else
        Ufield = sim.sol.u[i]
        x = sim.params.x
        Φr = @views Ufield[:, 3]
        q = sim.params.q
        xp1, vp1, _ = ParticleMotion.oscillator(sim.sol.t[i], sim.params.x0,
                                                sim.params.A, sim.params.ω)

        interpolator_Φr = QuadraticInterpolation(Φr, x)
        Φrz1 = interpolator_Φr(xp1)

        energy += -q * Φrz1
        return energy
    end
end

function FieldEnergy(sim)
    nt = length(sim.sol.t)
    energy = zeros(nt)

    for i in 1:nt
        energy[i] = Hamiltonian(i, sim)
    end

    return energy # ./ energy[1]
end

function MomentumVector_1()
    p = (Πr + Πs1)
    T0x = p * (Ψr + Ψs1)
end

function MomentumVector_2()
    p = (Πr + Πs1 + Πs2)
    T0x = p * (Ψr + Ψs1 + Ψs2)
end

function FieldMomentum_1()
    newton_cotes(MomentumVector_1)
end

end # end of module
