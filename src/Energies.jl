module Energies

using ..ScalarField
using ..ParticleMotion
using DataInterpolations

function newton_cotes(y::AbstractVector{T}, h::T; order::Int=1) where {T}
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

function regular_hamiltonian_density(i, sim)
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

    Hr = @. 0.5 * (Πr * Πr + Ψr * Ψr)
    return Hr
end

function interaction_hamiltonian_density(i, sim)
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

    # Hi = @. -Ψr * (Ψs1 + Ψs2)
    Hi = @. Ψr * (Ψs1 + Ψs2) + Πr * (Πs1 + Πs2)
    return Hi
end

function singular_hamiltonian_density(i, sim)
    Ufield = sim.sol.u[i].x[1]
    Upart = sim.sol.u[i].x[2]

    x = sim.params.x
    q1 = sim.params.q1
    q2 = sim.params.q2

    xp1 = Upart[2]
    xp2 = Upart[5]
    vp1 = Upart[3]
    vp2 = Upart[6]

    Πs1 = ScalarField.Πs.(x, q1, xp1, vp1)
    Ψs1 = ScalarField.Ψs.(x, q1, xp1, vp1)

    Πs2 = ScalarField.Πs.(x, q2, xp2, vp2)
    Ψs2 = ScalarField.Ψs.(x, q2, xp2, vp2)

    Hs = @. 0.5(Πs1 + Πs2)^2 + 0.5(Ψs1 + Ψs2)^2
    return Hs
end

function hamiltonian_density(i, sim)
    Hr = regular_hamiltonian_density(i, sim)
    Hi = interaction_hamiltonian_density(i, sim)
    # Hs = singular_hamiltonian_density(i, sim)
    H = Hr + Hi #+ Hs
    return H
end

function hamiltonian_density_1(i, sim)
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

function regular_hamiltonian(i, sim)
    H = regular_hamiltonian_density(i, sim)
    energy = newton_cotes(H, sim.params.dx; order=1)

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

    interpolator_Φr = PCHIPInterpolation(Φr, x)
    Φrz1 = interpolator_Φr(xp1)
    Φrz2 = interpolator_Φr(xp2)
    tmp = -q1 * Φrz1 - q2 * Φrz2 - q1 * Φs2_at_1 - q2 * Φs1_at_2
    # @show tmp
    # energy -= -q1 * Φrz1 - q2 * Φrz2 - q1 * Φs2_at_1 - q2 * Φs1_at_2
    return energy
end

function interaction_hamiltonian(i, sim)
    H = interaction_hamiltonian_density(i, sim)
    energy = newton_cotes(H, sim.params.dx; order=1)
    return energy
end

function singular_hamiltonian(i, sim)
    H = singular_hamiltonian_density(i, sim)
    energy = newton_cotes(H, sim.params.dx; order=1)
    Upart = sim.sol.u[i].x[2]

    q1 = sim.params.q1
    q2 = sim.params.q2

    xp1 = Upart[2]
    xp2 = Upart[5]
    vp1 = Upart[3]
    vp2 = Upart[6]

    Φs1_at_2 = ScalarField.Φs.(xp2, q1, xp1, vp1)
    Φs2_at_1 = ScalarField.Φs.(xp1, q2, xp2, vp2)
    # energy += -q1 * Φs2_at_1 - q2 * Φs1_at_2
    return energy
end

function hamiltonian(i, sim)
    Hr = regular_hamiltonian(i, sim)
    Hi = interaction_hamiltonian(i, sim)
    # Hs = singular_hamiltonian(i, sim)
    return Hr + Hi # + Hs
end

function field_energy(sim)
    nt = length(sim.sol.t)
    energy = zeros(nt)

    for i in 1:nt
        energy[i] = hamiltonian(i, sim)
    end

    return energy # ./ energy[1]
end

function regular_momentum_density(i, sim)
    Ufield = sim.sol.u[i].x[1]
    Πr = @views Ufield[:, 1]
    Ψr = @views Ufield[:, 2]

    Tr = -Πr .* Ψr
    return Tr
end

function interaction_momentum_density(i, sim)
    Ufield = sim.sol.u[i].x[1]
    Upart = sim.sol.u[i].x[2]
    x = sim.params.x
    q1 = sim.params.q1
    q2 = sim.params.q2
    xp1 = Upart[2]
    xp2 = Upart[5]
    vp1 = Upart[3]
    vp2 = Upart[6]
    # Φs1 = ScalarField.Φs.(x, q1, xp1, vp1)
    Πs1 = ScalarField.Πs.(x, q1, xp1, vp1)
    Ψs1 = ScalarField.Ψs.(x, q1, xp1, vp1)

    # Φs2 = ScalarField.Φs.(x, q2, xp2, vp2)
    Πs2 = ScalarField.Πs.(x, q2, xp2, vp2)
    Ψs2 = ScalarField.Ψs.(x, q2, xp2, vp2)
    Πr = @views Ufield[:, 1]
    Ψr = @views Ufield[:, 2]

    Ti = @. -Πs1 * Ψr - Πs2 * Ψr - Πr * Ψs1 - Πr * Ψs2
    return Ti
end

function singular_momentum_density(i, sim)
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
    Ts = @. -Πs1 * Ψs1 - Πs2 * Ψs1 - Πs1 * Ψs2 - Πs2 * Ψs2
    return Ts
end

function momentum_density(i, sim)
    Tr = regular_momentum_density(i, sim)
    Ti = interaction_momentum_density(i, sim)
    Ts = singular_momentum_density(i, sim)
    return Tr + Ti + Ts
end

function field_momentum(i, sim)
    T = momentum_density(i, sim)
    P = newton_cotes(T, sim.params.dx; order=1)
    return P
end

function half_field_momentum(i, sim)
    T = momentum_density(i, sim)
    N = length(T)
    t = @view T[div(N, 2):end]
    P = newton_cotes(t, sim.params.dx; order=1)
    return P
end

function half_field_regular_momentum(i, sim)
    T = regular_momentum_density(i, sim)
    N = length(T)
    t = @view T[div(N, 2):end]
    P = newton_cotes(t, sim.params.dx; order=1)
    return P
end

function half_field_interaction_momentum(i, sim)
    T = interaction_momentum_density(i, sim)
    N = length(T)
    t = @view T[div(N, 2):end]
    P = newton_cotes(t, sim.params.dx; order=1)
    return P
end

function half_field_singular_momentum(i, sim)
    T = singular_momentum_density(i, sim)
    N = length(T)
    t = @view T[div(N, 2):end]
    P = newton_cotes(t, sim.params.dx; order=1)
    return P
end

function total_hamiltonian_density(i, sim)
    Ufield = sim.sol.u[i].x[1]
    Upart = sim.sol.u[i].x[2]
    x = sim.params.x

    Πr = @views Ufield[:, 1]
    Ψr = @views Ufield[:, 2]
    Φr = @views Ufield[:, 3]

    q1 = sim.params.q1
    q2 = sim.params.q2

    mp1 = Upart[1]
    xp1 = Upart[2]
    vp1 = Upart[3]

    mp2 = Upart[4]
    xp2 = Upart[5]
    vp2 = Upart[6]

    Φs1(x) = ScalarField.Φs(x, q1, xp1, vp1)
    Φs2(x) = ScalarField.Φs(x, q2, xp2, vp2)
    Φs(x) = Φs1(x) + Φs2(x)

    Πs1(x) = ScalarField.Πs(x, q1, xp1, vp1)
    Πs2(x) = ScalarField.Πs(x, q2, xp2, vp2)
    Πs(x) = Πs1(x) + Πs2(x)

    Ψs1(x) = ScalarField.Ψs(x, q1, xp1, vp1)
    Ψs2(x) = ScalarField.Ψs(x, q2, xp2, vp2)
    Ψs(x) = Ψs1(x) + Ψs2(x)

    Πs1v1(x) = ScalarField.∂vΠ(x, q1, xp1, vp1)
    Πs2v2(x) = ScalarField.∂vΠ(x, q2, xp2, vp2)

    Ψs1v1(x) = ScalarField.∂vΨ(x, q1, xp1, vp1)
    Ψs2v2(x) = ScalarField.∂vΨ(x, q2, xp2, vp2)

    # Φs1_at_2 = ScalarField.Φs.(xp2, q1, xp1, vp1)
    # Φs2_at_1 = ScalarField.Φs.(xp1, q2, xp2, vp2)

    # interpolator_Φr = PCHIPInterpolation(Φr, x)
    # Φrz1 = interpolator_Φr(xp1)
    # Φrz2 = interpolator_Φr(xp2)

    Hfr = @. 0.5 * (Πr * Πr + Ψr * Ψr)
    Hfi = @. Πr * Πs.(x) + Ψr * Ψs.(x)
    Hfs = @. 0.5 * (Πs.(x) * Πs.(x) + Ψs.(x) * Ψs.(x))

    Hp1 = @. (Πr + Πs.(x)) * vp1 * Πs1v1.(x) - (Ψr + Ψs.(x)) * vp1 * Ψs1v1.(x)
    Hp2 = @. (Πr + Πs.(x)) * vp2 * Πs2v2.(x) - (Ψr + Ψs.(x)) * vp2 * Ψs2v2.(x)
    return Hfr, Hfi, Hfs, Hp1, Hp2
end

### Particle
###

@inline function KineticEnergy(m, v)
    γ = one(v) / sqrt(one(v) - v * v)
    return m * (γ - one(γ))
end

@inline function PotentialEnergy(xp, xs, Φ)
    interpolator_Φ = PCHIPInterpolation(Φ, xs)
    return interpolator_Φ(xp)
end

@inline function TotalEnergy(m, v, xp, xs, Φ)
    return m + KineticEnergy(m, v) + PotentialEnergy(xp, xs, Φ)
end

end # end of module
