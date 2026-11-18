module Diagnostics

using ..Equations
using ..ODE
using ..ScalarField
using StaticArrays: SVector

export semidiscrete_energy_terms, particle_observables, runtime_health,
       remainder_energy_balance

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

@inline function _cumulative_trapezoid(t, y)
    out = zeros(promote_type(eltype(t), eltype(y)), length(t))
    for i in 2:length(t)
        out[i] = out[i - 1] + 0.5 * (t[i] - t[i - 1]) * (y[i - 1] + y[i])
    end
    return out
end

"""Evaluate the actual periodic semidiscrete remainder-energy balance over a solution."""
function remainder_energy_balance(sim)
    t = sim.sol.t
    params = sim.sol.prob.p
    h = params.h
    c = params.equation.velocity
    nt = length(t)
    energy = zeros(nt)
    source_power = zeros(nt)
    dissipation = zeros(nt)
    reconstruction = zeros(nt)

    for j in eachindex(t)
        Q = sim.sol.u[j].x[1]
        particle = sim.sol.u[j].x[2]
        m, z, v = particle[1], particle[2], particle[3]
        _, Πr, Ψr = ODE.regular_fields_at_particle(Q, z, params)
        feedback = hasproperty(params, :particle_feedback) ? params.particle_feedback : true
        rates = ODE.particle_rates(m, z, v, params.q, Πr, Ψr,
                                   params.dxPotential;
                                   particle_feedback=feedback)
        source = zeros(eltype(Q), size(Q))
        if !hasproperty(params, :field_source) || params.field_source
            for i in 1:params.N
                ODE.add_singular_sources!(source, i, params.x[i], h, z, v,
                                          params.q, rates.dv)
            end
        end
        terms = semidiscrete_energy_terms(Q, source, params)
        energy[j] = 0.5h * sum(@view(Q[:, 1]) .^ 2 ./ c^2 .+
                               @view(Q[:, 2]) .^ 2)
        source_power[j] = terms.source_power
        dissipation[j] = terms.dissipation
        reconstruction[j] = terms.reconstruction_residual
    end

    integrated_source = _cumulative_trapezoid(t, source_power)
    integrated_dissipation = _cumulative_trapezoid(t, dissipation)
    integrated_reconstruction = _cumulative_trapezoid(t, reconstruction)
    energy_change = energy .- first(energy)
    predicted_change = integrated_source .- integrated_dissipation .+
                       integrated_reconstruction
    residual = energy_change .- predicted_change
    return (; t, energy, energy_change, source_power, dissipation, reconstruction,
            integrated_source, integrated_dissipation, integrated_reconstruction,
            predicted_change, residual)
end

end
