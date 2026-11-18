module Convergence

using ..Run
using ..InitialData
using ..PlottingUtils: mytheme_aps
using CairoMakie
using DataInterpolations

function convergence(L, N1, tf=1.0, cfl=0.5, sf=true, muscl=true)
    x1, sol1 = ScalarWaveFVM.Run.run(L, N1, tf, cfl, sf, muscl)
    x2, sol2 = ScalarWaveFVM.Run.run(L, 2N1, tf, cfl, sf, muscl)
    x4, sol4 = ScalarWaveFVM.Run.run(L, 4N1, tf, cfl, sf, muscl)
    return x1, x2, x4, sol1, sol2, sol4
end

function run_coupled_sims(L, N, tf, cfl, sf)
    sim1 = Run.coupled_system(L, N, tf, cfl, sf)
    sim2 = Run.coupled_system(L, 2N, tf, cfl, sf)
    sim4 = Run.coupled_system(L, 4N, tf, cfl, sf)
    return sim1, sim2, sim4
end

function run_forced_sims(Nosc, dx0, cfl, vmax, q; sf=true, muscl=true)
    sim1 = Run.forced_motion(Nosc, dx0, cfl, vmax, q, sf, muscl)
    sim2 = Run.forced_motion(Nosc, dx0 / 2, cfl, vmax, q, sf, muscl)
    sim4 = Run.forced_motion(Nosc, dx0 / 4, cfl, vmax, q, sf, muscl)
    return sim1, sim2, sim4
end

function mymean(arr)
    arr2 = zeros(eltype(arr), div(length(arr), 2))
    for i in eachindex(arr2)
        j = 2i - 1
        arr2[i] = 0.5(arr[j] + arr[j + 1])
    end
    return arr2
end

function total_variation(q::AbstractArray)
    TV = zero(eltype(q))
    for i in eachindex(q)
        TV += abs(q[i] - q[i - 1])
    end
    return TV
end

function convergence_order(p, x1, x2, sol1, sol2, true_sol)
    order = zeros(length(sol1.t))
    for i in 1:(length(sol1.t))
        t1 = sol1.t[i]
        t2 = sol2.t[2i - 1]
        t1 == t2 || println(t1, "   ", t2)
        t = t1
        u1 = true_sol.(t, x1)
        u2 = true_sol.(t, x2)
        tmp1 = sum(abs.(sol1[:, 1, i] .- u1) .^ p)^(1 / p)
        tmp2 = sum(abs.(mymean(sol2[:, 1, 2i - 1] .- u2)) .^ p)^(1 / p)
        order[i] = log2(tmp1 / tmp2)
    end
    return order
end

function self_convergence_order(p, x, sol1, sol2, sol4)
    order = zeros(length(sol1.t))
    for i in 1:(length(sol1.t))
        tmp1 = sum(abs.(sol1[:, 1, i] .- mymean(sol2[:, 1, 2i - 1])) .^ p)^(1 / p)
        tmp2 = sum(abs.(mymean(sol2[:, 1, 2i - 1]) .- mymean(mymean(sol4[:, 1, 4i - 3]))) .^
                   p)^(1 /
                       p)
        order[i] = log2(tmp1 / tmp2)
    end
    return order
end

function particle_energy(sim)
    ts = sim.sol.t
    nt = length(ts)
    x = sim.params.x
    q = sim.params.q
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    interpolation_method = sim.sol.prob.p.interpolation_method
    V(xx) = InitialData.Harmonic(xx, k)

    energy = zeros(nt)
    for i in eachindex(ts)
        z = sim.sol.u[i].x[2][2]
        v = sim.sol.u[i].x[2][3]
        γ = inv(sqrt(1.0 - v * v))
        Φr = @view sim.sol.u[i].x[1][:, 3]
        interpolator_Φ = interpolation_method(Φr, x;
                                              extrapolation=ExtrapolationType.Extension)
        Φrz = interpolator_Φ(z)
        energy[i] = (m0 - q * Φrz) * (γ - 1.0) + V(z)
    end
    return energy
end

function no_radiation_particle_energy(sim)
    ts = sim.sol.t
    nt = length(ts)
    m0 = hasproperty(sim.params, :m0) ? sim.params.m0 : sim.params.m
    k = sim.params.k
    V(xx) = InitialData.Harmonic(xx, k)
    sol = Run.solve_non_radiative_relativistic_particle(sim.params)

    energy = zeros(nt)
    for i in eachindex(ts)
        z = sol.u[i][2]
        v = sol.u[i][3]
        γ = inv(sqrt(1.0 - v * v))
        energy[i] = m0 * (γ - 1.0) + V(z)
    end
    return energy
end

function particle_energy_convergence_to_no_radiation(sim1, sim2, sim4)
    E1 = particle_energy(sim1)
    E2 = particle_energy(sim2)
    E4 = particle_energy(sim4)

    Enr1 = no_radiation_particle_energy(sim1)
    Enr2 = no_radiation_particle_energy(sim2)
    Enr4 = no_radiation_particle_energy(sim4)

    nt = length(sim1.sol.t)

    error1 = zeros(nt)
    error2 = zeros(nt)
    error4 = zeros(nt)
    order12 = zeros(nt)
    order24 = zeros(nt)
    order = zeros(nt)

    for i in 1:nt
        error1[i] = abs(E1[i] - Enr4[4i - 3])#Enr1)
        error2[i] = abs(E2[2i - 1] - Enr4[4i - 3])#Enr2[2i - 1])
        error4[i] = abs(E4[4i - 3] - Enr4[4i - 3])
        order12[i] = log2((error1[i] + eps(Float64)) / (error2[i] + eps(Float64)))
        order24[i] = log2((error2[i] + eps(Float64)) / (error4[i] + eps(Float64)))
        order[i] = 0.5 * (order12[i] + order24[i])
    end

    return (; t=sim1.sol.t, error1, error2, error4, order12, order24, order)
end

function plot_particle_energy_convergence_to_no_radiation(sim1, sim2, sim4)
    conv = particle_energy_convergence_to_no_radiation(sim1, sim2, sim4)

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              title=L"\textrm{Particle-energy convergence to no-radiation}")

    lines!(ax, conv.t[2:end], conv.order12[2:end];
           label=L"\log_2(e_1/e_2)", linewidth=1.5)
    lines!(ax, conv.t[2:end], conv.order24[2:end];
           label=L"\log_2(e_2/e_4)", linewidth=1.5, linestyle=:dash)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)
    return fig
end

function plot_particle_energy_rescaled_errors_to_no_radiation(sim1, sim2, sim4)
    conv = particle_energy_convergence_to_no_radiation(sim1, sim2, sim4)
    eps0 = eps(Float64)

    e1 = max.(conv.error1, eps0)
    e2 = max.(conv.error2, eps0)
    e4 = max.(conv.error4, eps0)

    CairoMakie.activate!()
    mytheme = mytheme_aps()
    set_theme!(mytheme)

    fig = Figure(; size=(253, 200), figure_padding=(1, 5, 5, 1))
    ax = Axis(fig[1, 1];
              xlabel=L"t",
              ylabel=L"\textrm{rescaled error}",
              yscale=log10,
              title=L"\textrm{Rescaled particle-energy errors}")
    ylims!(ax, 1e-5, 1e-2)
    lines!(ax, conv.t, e1; label=L"e_1", linewidth=1.5)
    lines!(ax, conv.t, 2.0 .* e2; label=L"2e_2", linewidth=1.5, linestyle=:dash)
    lines!(ax, conv.t, 4.0 .* e4; label=L"4e_4", linewidth=1.5, linestyle=:dot)

    Legend(fig[2, 1], ax; orientation=:horizontal, tellwidth=false, tellheight=true)
    rowsize!(fig.layout, 1, Auto())
    rowsize!(fig.layout, 2, Relative(0.12))
    colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 0)
    return fig
end

end # end of module
