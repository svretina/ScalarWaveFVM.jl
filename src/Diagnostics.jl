module Diagnostics

using ..Equations
using ..ODE
using ..ScalarField
using StaticArrays: SVector

export semidiscrete_energy_terms, particle_observables, runtime_health

@inline _energy_dot(q, p, c) = q[1] * p[1] / c^2 + q[2] * p[2]

"""
    semidiscrete_energy_terms(Q, source, params)

Return the exact periodic semidiscrete energy decomposition for the MUSCL
operator used by `ODE.muscl_optimized!`:

`rate = source_power - dissipation + reconstruction_residual`.

The energy is `h/2 Σ(Π²/c² + Ψ²)`. Godunov and Rusanov both have wave-speed
dissipation `c/2 ΔQᵀHΔQ`; the residual measures the exact change caused by
feeding reconstructed rather than cell-average states to the numerical flux.
"""
function semidiscrete_energy_terms(Q, source, params)
    N, h = params.N, params.h
    equation = params.equation
    c = equation.velocity
    flux = params.numerical_flux
    L = Equations.left_eigenmatrix(equation)
    R = Equations.right_eigenmatrix(equation)
    limiter = params.slope_limiter

    source_power = zero(eltype(Q))
    dissipation = zero(eltype(Q))
    reconstruction_residual = zero(eltype(Q))

    for i in 1:N
        ip1 = mod(i, N) + 1
        Qi = SVector(Q[i, 1], Q[i, 2])
        Qip1 = SVector(Q[ip1, 1], Q[ip1, 2])
        _, qL = ODE.reconstruct_cell(Q, i, N, L, R, limiter)
        qR, _ = ODE.reconstruct_cell(Q, i + 1, N, L, R, limiter)
        actual_flux = flux(qL, qR, equation)
        ΔQ = Qip1 - Qi
        cell_flux = Equations.principal_part(equation) * (Qi + Qip1) / 2 - c * ΔQ / 2

        dissipation += 0.5 * c * _energy_dot(ΔQ, ΔQ, c)
        reconstruction_residual += _energy_dot(ΔQ, actual_flux - cell_flux, c)
        source_power += h * _energy_dot(Qi,
                                        SVector(source[i, 1], source[i, 2]), c)
    end

    rate = source_power - dissipation + reconstruction_residual
    return (; rate, source_power, dissipation, reconstruction_residual,
            boundary_flux=zero(rate))
end

function particle_observables(Q_field, particle, params; m0=params.m0, τ0=0.0)
    m, z, v, τ = particle
    Φr, Πr, Ψr = ODE.regular_fields_at_particle(Q_field, z, params)
    reconstructed_mass = ODE.dynamical_mass(m0, params.q, Φr, τ; τ0=τ0)
    γ = inv(sqrt(1 - v^2))
    return (; m, z, v, τ, γ, Φr, Πr, Ψr, reconstructed_mass,
            mass_constraint=m - reconstructed_mass,
            scalar_force_combination=Ψr + v * Πr)
end

function runtime_health(states)
    masses = (u.x[2][1] for u in states)
    velocities = (u.x[2][3] for u in states)
    minimum_mass = minimum(masses)
    maximum_speed = maximum(abs, velocities)
    finite = all(u -> all(isfinite, u.x[1]) && all(isfinite, u.x[2]), states)
    return (; minimum_mass, maximum_speed, finite,
            valid=finite && minimum_mass > 0 && maximum_speed < 1)
end

end
