using ScalarWaveFVM
using CairoMakie
using DelimitedFiles
using LaTeXStrings

const RESOLUTIONS = (0.0500, 0.0250, 0.0125)
const N_OSCILLATIONS = 4.75
const CFL = 0.4
const V_MAX = 0.30
const CHARGE = 0.35

function particle_history(sim)
    particle = [u.x[2] for u in sim.sol.u]
    return (; t=sim.sol.t,
            mass=getindex.(particle, 1),
            position=getindex.(particle, 2),
            velocity=getindex.(particle, 3),
            proper_time=getindex.(particle, 4))
end

function turning_point_amplitudes(position)
    return [abs(position[i]) for i in 2:(length(position) - 1)
            if abs(position[i]) >= abs(position[i - 1]) &&
               abs(position[i]) > abs(position[i + 1])]
end

function interpolate_at(target_time, time, values)
    result = similar(target_time)
    for (i, target) in pairs(target_time)
        right = searchsortedfirst(time, target)
        if right <= 1
            result[i] = first(values)
        elseif right > length(time)
            result[i] = last(values)
        else
            left = right - 1
            weight = (target - time[left]) / (time[right] - time[left])
            result[i] = (1 - weight) * values[left] + weight * values[right]
        end
    end
    return result
end

sims = [ScalarWaveFVM.Run.particle_in_potential_relativistic(
            N_OSCILLATIONS, dx, CFL, V_MAX, CHARGE;
            amplitude=0.10,
            initial_position=0.10,
            initial_velocity=0.0,
            rest_mass=1.0,
        ) for dx in RESOLUTIONS]
histories = particle_history.(sims)

ScalarWaveFVM.enable_plotting!()
output_dir = joinpath(pkgdir(ScalarWaveFVM), "figures")
mkpath(output_dir)
stem = joinpath(output_dir, "damped_trajectory_resolution_comparison")
labels = [L"\Delta x=0.05", L"\Delta x=0.025", L"\Delta x=0.0125"]
figure = ScalarWaveFVM.ParticleHistoryPlots.plot_particle_trajectory_resolutions(
    sims, labels; savepath=stem * ".png")
save(stem * ".pdf", figure)

finest = last(histories)
summary = zeros(length(sims), 9)
for (i, (dx, sim, history)) in enumerate(zip(RESOLUTIONS, sims, histories))
    peaks = turning_point_amplitudes(history.position)
    finest_position = interpolate_at(history.t, finest.t, finest.position)
    maximum_difference = maximum(abs, history.position .- finest_position)
    summary[i, :] .= (dx, sim.sol.prob.p.h, sim.sol.prob.p.dt,
                      minimum(history.mass), maximum(abs, history.velocity),
                      first(peaks), last(peaks),
                      1 - last(peaks) / first(peaks), maximum_difference)

    history_data = hcat(history.t, history.mass, history.position,
                        history.velocity, history.proper_time)
    filename = stem * "_dx_$(replace(string(dx), '.' => 'p')).csv"
    open(filename, "w") do io
        println(io, "t,mass,position,velocity,proper_time")
        writedlm(io, history_data, ',')
    end
end

open(stem * "_summary.csv", "w") do io
    println(io, "requested_dx,actual_dx,dt,min_mass,max_abs_velocity," *
                "first_peak,last_peak,fractional_peak_decay," *
                "max_position_difference_vs_finest")
    writedlm(io, summary, ',')
end

println("SAVED=$(stem).png")
println("SUMMARY_ROWS=")
display(summary)
