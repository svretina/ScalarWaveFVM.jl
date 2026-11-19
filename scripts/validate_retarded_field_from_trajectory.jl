using ScalarWaveFVM
using CairoMakie
using DelimitedFiles
using LaTeXStrings

include(joinpath(pkgdir(ScalarWaveFVM), "validation", "RetardedFieldReference.jl"))
using .RetardedFieldReference

const RESOLUTIONS = (0.0500, 0.0250, 0.0125)
const VALIDATION_TIMES = (2.0, 4.0, 6.0, 8.0, 9.0)
const N_OSCILLATIONS = 4.75
const CFL = 0.4
const V_MAX = 0.30
const CHARGE = 0.35
const AMPLITUDE = 0.10
const REST_MASS = 1.0
const PLOT_ONLY = "--plot-only" in ARGS
const OUTPUT_DIR = joinpath(pkgdir(ScalarWaveFVM), "figures")
const STEM = joinpath(OUTPUT_DIR, "retarded_field_from_coupled_trajectory")

function trajectory_history(sim)
    particle = [state.x[2] for state in sim.sol.u]
    return TrajectoryHistory(collect(sim.sol.t),
                             getindex.(particle, 2),
                             getindex.(particle, 3),
                             getindex.(particle, 4))
end

function state_at_saved_time(sim, time)
    index = findfirst(t -> isapprox(t, time; atol=32eps(time)), sim.sol.t)
    index === nothing && error("time $time is not present in the saved solution")
    return sim.sol.u[index]
end

function compare_field(sim, history, time)
    params = sim.sol.prob.p
    h, x, L = params.h, params.x, params.L
    numerical = state_at_saved_time(sim, time).x[1]
    initial = trajectory_state(history, first(history.t))

    boundary_safe = @. L - abs(x) > time + 4h
    front_safe = @. abs(abs(x - initial.z) - time) > 4h
    mask = boundary_safe .& front_safe
    indices = findall(mask)
    reference = [cell_average_exact_remainder_1d(history, time, x[i], h, CHARGE)
                 for i in indices]
    reference_Π = getproperty.(reference, :Π)
    reference_Ψ = getproperty.(reference, :Ψ)
    reference_Φ = getproperty.(reference, :Φ)
    error_Π = numerical[indices, 1] .- reference_Π
    error_Ψ = numerical[indices, 2] .- reference_Ψ
    error_Φ = numerical[indices, 3] .- reference_Φ

    metrics = (L1_Π=h * sum(abs, error_Π),
               L1_Ψ=h * sum(abs, error_Ψ),
               L1_Φ=h * sum(abs, error_Φ),
               Linf_Π=maximum(abs, error_Π),
               Linf_Ψ=maximum(abs, error_Ψ),
               Linf_Φ=maximum(abs, error_Φ),
               cells=length(indices))
    profile = (; x=x[indices],
               numerical_Π=numerical[indices, 1], reference_Π, error_Π,
               numerical_Ψ=numerical[indices, 2], reference_Ψ, error_Ψ,
               numerical_Φ=numerical[indices, 3], reference_Φ, error_Φ)
    return metrics, profile
end

mkpath(OUTPUT_DIR)
if PLOT_ONLY
    data = readdlm(STEM * "_fine_profile.csv", ',', Float64; skipstart=1)
    final_profile = (x=data[:, 1],
                     numerical_Π=data[:, 2], reference_Π=data[:, 3],
                     error_Π=data[:, 4],
                     numerical_Ψ=data[:, 5], reference_Ψ=data[:, 6],
                     error_Ψ=data[:, 7],
                     numerical_Φ=data[:, 8], reference_Φ=data[:, 9],
                     error_Φ=data[:, 10])
    summary = readdlm(STEM * "_summary.csv", ',', Float64; skipstart=1)
else
    rows = Vector{Vector{Float64}}()
    final_profile = nothing
    for (resolution_index, requested_dx) in enumerate(RESOLUTIONS)
        sim = ScalarWaveFVM.Run.particle_in_potential_relativistic(
            N_OSCILLATIONS, requested_dx, CFL, V_MAX, CHARGE;
            amplitude=AMPLITUDE,
            initial_position=AMPLITUDE,
            initial_velocity=0.0,
            rest_mass=REST_MASS,
        )
        history = trajectory_history(sim)
        for time in VALIDATION_TIMES
            metrics, profile = compare_field(sim, history, time)
            push!(rows, [requested_dx, sim.sol.prob.p.h, time, metrics.cells,
                         metrics.L1_Π, metrics.L1_Ψ, metrics.L1_Φ,
                         metrics.Linf_Π, metrics.Linf_Ψ, metrics.Linf_Φ])
            if resolution_index == length(RESOLUTIONS) &&
               time == last(VALIDATION_TIMES)
                global final_profile = profile
            end
        end
        sim = nothing
        history = nothing
        GC.gc()
    end
    summary = reduce(vcat, permutedims.(rows))
    open(STEM * "_summary.csv", "w") do io
        println(io, "requested_dx,actual_dx,time,compared_cells," *
                    "L1_Pi,L1_Psi,L1_Phi,Linf_Pi,Linf_Psi,Linf_Phi")
        writedlm(io, summary, ',')
    end

    profile_data = hcat(final_profile.x,
                        final_profile.numerical_Π, final_profile.reference_Π,
                        final_profile.error_Π,
                        final_profile.numerical_Ψ, final_profile.reference_Ψ,
                        final_profile.error_Ψ,
                        final_profile.numerical_Φ, final_profile.reference_Φ,
                        final_profile.error_Φ)
    open(STEM * "_fine_profile.csv", "w") do io
        println(io, "x,numerical_Pi,reference_Pi,error_Pi," *
                    "numerical_Psi,reference_Psi,error_Psi," *
                    "numerical_Phi,reference_Phi,error_Phi")
        writedlm(io, profile_data, ',')
    end
end

ScalarWaveFVM.enable_plotting!()
set_theme!(ScalarWaveFVM.PlottingUtils.mytheme_aps())
fig = Figure(; size=(510, 520), figure_padding=(1, 5, 5, 1))
fields = ((final_profile.numerical_Π, final_profile.reference_Π,
           final_profile.error_Π, L"\Pi_r"),
          (final_profile.numerical_Ψ, final_profile.reference_Ψ,
           final_profile.error_Ψ, L"\Psi_r"),
          (final_profile.numerical_Φ, final_profile.reference_Φ,
           final_profile.error_Φ, L"\Phi_r"))
for (row, (numerical, reference, error, label)) in enumerate(fields)
    value_axis = Axis(fig[row, 1]; ylabel=label,
                      title=row == 1 ? L"\textrm{Numerical and Retarded Reference}" : "")
    error_axis = Axis(fig[row, 2]; ylabel=L"U_{r,h}-\overline{U}_r^{\mathrm{ref}}",
                      title=row == 1 ? L"\textrm{Cell-Average Error}" : "")
    lines!(value_axis, final_profile.x, numerical; label=L"\textrm{numerical}")
    lines!(value_axis, final_profile.x, reference;
           linestyle=:dash, label=L"\textrm{reference}")
    lines!(error_axis, final_profile.x, error)
    row == 1 && axislegend(value_axis; position=:rb, labelsize=8, rowgap=0.1)
    if row < 3
        hidexdecorations!(value_axis; grid=false)
        hidexdecorations!(error_axis; grid=false)
    else
        value_axis.xlabel = L"x"
        error_axis.xlabel = L"x"
    end
end
rowgap!(fig.layout, 4)
colgap!(fig.layout, 8)
save(STEM * ".png", fig)
save(STEM * ".pdf", fig)

println("SAVED=$(STEM).png")
println("SUMMARY_ROWS=")
display(summary)
