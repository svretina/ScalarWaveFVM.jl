module MyPlots

using Plots
using LaTeXStrings
using ..ParticleMotion
using ..Run

Plots.PGFPlotsXBackend()
default(; legend=:topright,              # Legend position
        framestyle=:box,              # Boxed frame style for axes
        linewidth=1.5,                # Line width for plot lines
        grid=false,                   # Disable grid lines (as seen in the PDF)
        label="",                     # No labels for the lines by default
        tickfontsize=20,              # Font size for tick labels
        guidefontsize=20,             # Font size for axis labels
        legendfontsize=12,            # Font size for legend text
        foreground_color_border=:black, # Black axis lines
        size=(800, 600),              # Size of the plot
        fontfamily="Computer Modern",  # Use Computer Modern font
        dpi=200,
        thickness_scaling=2)

function mymean(arr)
    arr2 = zeros(eltype(arr), div(length(arr), 2))
    for i in eachindex(arr2)
        j = 2i - 1
        arr2[i] = 0.5(arr[j] + arr[j + 1])
    end
    return arr2
end

function plot_errors(i, true_sol, scale=2)
    t1 = sol1.t[i]
    t2 = sol2.t[2i - 1]
    @assert t1 == t2
    t = t1
    p = scatter(x1, sol1[:, 1, i] .- true_sol.(t, x1))
    scatter!(p, x1, (mymean(sol2[:, 1, 2i - 1]) .- true_sol.(t, x1)) .* scale)
    title!(p, "t=$t")
    return p
end

function plot_self_errors(ti, scale=2)
    plot1 = scatter(x1,
                    (mymean(sol2[:, 1, 2ti - 1]) .-
                     mymean(mymean(sol4[:, 1, 4ti - 3]))) .* scale)
    scatter!(plot1, x1, sol1[:, 1, ti] .- sol2[1:2:end, 1, 2ti - 1])
    return plot1
end

function plot_field_resolutions(i, sf=true)
    sim1 = Run.coupled_system(1000, 2000, 800, 0.25, sf)
    sim2 = Run.coupled_system(1000, 4000, 800, 0.25, sf)
    sim4 = Run.coupled_system(1000, 8000, 800, 0.25, sf)

    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1
    p = scatter(sim1.params.x, sim1.sol.u[i].x[1][:, 1]; label="low res")
    scatter!(p, sim1.params.x, mymean(sim2.sol.u[2i - 1].x[1][:, 1]); label="mid res")
    scatter!(p, sim1.params.x, mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 1]));
             label="high res")
    xaxis!(p, "x")
    yaxis!(p, "Π")
    title!(p, "t[$(i)]=$(round(t, digits=3))")
    savefig(p, "pifield_i=$(i).pdf")
    return p
end

function plot_errors_coupled(i, sf=true)
    sim1 = Run.coupled_system(1000, 2000, 800, 0.25, sf)
    sim2 = Run.coupled_system(1000, 4000, 800, 0.25, sf)
    sim4 = Run.coupled_system(1000, 8000, 800, 0.25, sf)

    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1
    p = plot()
    scatter!(p, sim1.params.x,
             sim1.sol.u[i].x[1][:, 1] .- mymean(sim2.sol.u[2i - 1].x[1][:, 1]);
             label="low-mid")
    scatter!(p, sim1.params.x,
             mymean(sim2.sol.u[2i - 1].x[1][:, 1])
             .-
             mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 1]));
             label="mid-high")
    return p
end

function plot_field_resolutions(i, sim1, sim2, sim4)
    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1
    p = scatter(sim1.params.x, sim1.sol.u[i].x[1][:, 1]; label="low res")
    scatter!(p, sim1.params.x, mymean(sim2.sol.u[2i - 1].x[1][:, 1]); label="mid res")
    scatter!(p, sim1.params.x, mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 1]));
             label="high res")
    xaxis!(p, "x")
    yaxis!(p, "Π")
    title!(p, "t[$(i)]=$(round(t, digits=3))")
    savefig(p, "pifield_i=$(i).pdf")
    return p
end

function plot_field_errors(i, sim1, sim2, sim4, scale=1.0)
    t1 = sim1.sol.t[i]
    t2 = sim2.sol.t[2i - 1]
    t4 = sim4.sol.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1
    p = plot()
    scatter!(p, sim1.params.x,
             sim1.sol.u[i].x[1][:, 1] .- mymean(sim2.sol.u[2i - 1].x[1][:, 1]);
             label="low-mid")
    scatter!(p, sim1.params.x,
             scale .* (mymean(sim2.sol.u[2i - 1].x[1][:, 1]) .-
                       mymean(mymean(sim4.sol.u[4i - 3].x[1][:, 1])));
             label="mid-high")
    return p
end

# function make_gif_errors(scale=1.0)
#     nt = length(sol1.t)
#     p = plot()
#     ylims!(p, -0.1, 1.1)
#     anim = @animate for i in 1:10:nt
#         plot_errors(i, scale)
#     end
#     gif(anim, "scaled_errors_fps15.gif"; fps=15)
# end

# function make_gif(x, sol)
#     nt = length(sol.t)
#     p = plot()
#     ylims!(p, -0.1, 1.1)
#     anim = @animate for i in 1:10:nt
#         scatter(p, x, sol[:, 1, i])
#     end
#     gif(anim, "anim_fps15.gif"; fps=15)
# end

# function make_gif_with_forced_particle(x, sol, params)
#     nt = length(sol.t)

#     x10 = params.x1
#     x20 = params.x2
#     L = params.L
#     direction1 = params.direction1
#     direction2 = params.direction2
#     vmax = params.vmax
#     f0 = params.f0
#     A = params.A
#     anim = @animate for i in 1:10:nt
#         p = plot()
#         ylims!(p, -3.1, 5.5)
#         t = sol.t[i]
#         x1, v1, a1 = ParticleMotion.oscillator2(t, x10, vmax, A, direction1)
#         x2, v2, a2 = ParticleMotion.oscillator2(t, x20, vmax, A, direction2)
#         plot!(p, x, sol[:, 1, i]; label=false)
#         vline!(p, [x1]; label=false)
#         vline!(p, [x2]; label=false)
#         title!(p, "Π at t[$(i)]=$(round(t,digits=3))")
#     end
#     gif(anim, "anim_fps15.gif"; fps=15)
# end

# function make_gif_with_particle(xs, sol)
#     nt = length(sol.t)
#     ymins = zeros(nt)
#     ymaxs = zeros(nt)
#     for i in 1:nt
#         ymins[i] = min(sol.u[i].x[1][:, 1]...)
#         ymaxs[i] = max(sol.u[i].x[1][:, 1]...)
#     end
#     ymin = min(ymins...)
#     ymax = max(ymaxs...)
#     p = plot()
#     anim = @animate for i in 1:10:nt
#         p = plot()
#         # ylims!(p, 1.1ymin, 1.1ymax)
#         t = sol.t[i]
#         x1 = sol.u[i].x[2][2]

#         plot!(p, xs, sol.u[i].x[1][:, 1]; label=false)
#         vline!(p, [x1]; label=false)
#         title!(p, "Π at t[$(i)]=$(round(t,digits=3))")
#     end
#     gif(anim, "anim_fps15.gif"; fps=15)
# end

# function make_gif_with_particles(x, sol)
#     nt = length(sol.t)
#     ymins = zeros(nt)
#     ymaxs = zeros(nt)
#     for i in 1:nt
#         ymins[i] = min(sol.u[i].x[1][:, 1]...)
#         ymaxs[i] = max(sol.u[i].x[1][:, 1]...)
#     end
#     ymin = min(ymins...)
#     ymax = max(ymaxs...)
#     anim = @animate for i in 1:10:nt
#         p = plot()
#         # ylims!(p, 1.1ymin, 1.1ymax)
#         t = sol.t[i]
#         x1 = sol.u[i].x[2][2]
#         x2 = sol.u[i].x[2][5]

#         plot!(p, x, sol.u[i].x[1][:, 1]; label=false)
#         vline!(p, [x1]; label=false)
#         vline!(p, [x2]; label=false)
#         title!(p, "Π at t[$(i)]=$(round(t,digits=3))")
#     end
#     gif(anim, "anim_fps15.gif"; fps=15)
# end

function plot_position(sol, existing_fig=nothing; invert=false)
    if existing_fig === nothing
        # Create a new figure with 1 row and 2 columns if no figure is provided
        p = plot(; layout=(1, 2))
    else
        # Use the existing figure
        p = existing_fig
        # Get the layout dimensions
        rows, cols = size(p.layout)
        # Check if there's a second row available
        if rows < 2
            error("The existing figure must have at least 2 rows to plot in the second row")
        end
    end
    if existing_fig === nothing
        idx1, idx2 = 1, 2
    else
        idx1 = cols + 1      # Start of second row
        idx2 = cols + 2      # Second position in second row
        # Check if these indices exist in the layout
        if idx2 > rows * cols
            error("Not enough subplots in the existing figure's layout for second row plotting")
        end
    end
    nt = length(sol.t)
    x1 = zeros(nt)
    x2 = zeros(nt)
    for (i, ti) in enumerate(sol.t)
        x1[i] = sol.u[i].x[2][2]
        x2[i] = sol.u[i].x[2][5]
    end
    if invert
        plot!(p[idx1], x1, sol.t; title="particle 1", xlabel="x", ylabel="t", label=false)
        plot!(p[idx2], x2, sol.t; title="particle 2", xlabel="x", ylabel="", label=false)
    else
        plot!(p[idx1], sol.t, x1; title="particle 1", xlabel="t", ylabel="x", label=false)
        plot!(p[idx2], sol.t, x2; title="particle 2", xlabel="t", ylabel="", label=false)
    end
    return p
end

function plot_position1(sim; invert=false)
    pgfplotsx()
    p = plot()

    nt = length(sim.sol.t)
    x1 = zeros(nt)
    for i in 1:nt
        x1[i] = sim.sol.u[i].x[2][2]
    end
    if invert
        plot!(p, x1, sim.sol.t; title="particle 1", xlabel="x", ylabel="t", label=false)
    else
        plot!(p, sim.sol.t, x1; title="particle 1", xlabel="t", ylabel="x", label=false)
    end
    savefig(p, "particle_location.pdf")
    return p
end

function plot_position_compare_sf(sim1, sim2, invert=false)
    pgfplotsx()

    @assert all(sim1.sol.t .== sim2.sol.t)
    p = plot(; title=L"m=%$(sim1.m),
q=%$(sim1.q), v_{0}=%$(sim1.v0), x_{0}=%$(sim1.x0)")

    nt1 = length(sim1.sol.t)
    nt2 = length(sim2.sol.t)
    x1 = zeros(nt1)
    x2 = zeros(nt2)

    for i in 1:nt1
        x1[i] = sim1.sol.u[i].x[2][2]
    end
    for i in 1:nt2
        x2[i] = sim2.sol.u[i].x[2][2]
    end

    if sim2.params.sf
        xtmp = copy(x1)
        x1 .= x2
        x2 .= xtmp
    end
    if invert
        plot!(p, x1, sim1.sol.t; xlabel="x", ylabel="t", label="self force")
        plot!(p, x2, sim2.sol.t; label="no self force")
    else
        plot!(p, sim1.sol.t, x1; xlabel="t", ylabel="x", label="self force")
        plot!(p, sim2.sol.t, x2; label="no self force")
    end
    savefig(p, "comparison_position.pdf")
    return p
end

function plot_position_resolutions(sim1, sim2, sim4; invert=false)
    pgfplotsx()
    p = plot()

    nt = length(sim1.sol.t)
    x1 = zeros(nt)
    x2 = zeros(nt)
    x4 = zeros(nt)
    for i in 1:nt
        x1[i] = sim1.sol.u[i].x[2][2]
        x2[i] = sim2.sol.u[2i - 1].x[2][2]
        x4[i] = sim4.sol.u[4i - 3].x[2][2]
    end
    if sim1.params.sf
        label_low = "low + sf"
        label_mid = "mid + sf"
        label_high = "high + sf"
    else
        label_low = "low"
        label_mid = "mid"
        label_high = "high"
    end
    if invert
        plot!(p, x1, sim1.sol.t; label=label_low)
        plot!(p, x2, sim1.sol.t; label=label_mid)
        plot!(p, x4, sim1.sol.t; label=label_high)
    else
        plot!(p, sim1.sol.t, x1; label=label_low)
        plot!(p, sim1.sol.t, x2; label=label_mid)
        plot!(p, sim1.sol.t, x4; label=label_high)
    end
    savefig(p, "particle_location_resolutions.pdf")
    return p
end

function plot_position_errors(sim1, sim2, sim4, scale=2; invert=false)
    pgfplotsx()
    p = plot(; title="Errors")

    nt = length(sim1.sol.t)
    x1 = zeros(nt)
    x2 = zeros(nt)
    x4 = zeros(nt)
    for i in 1:nt
        x1[i] = sim1.sol.u[i].x[2][2]
        x2[i] = sim2.sol.u[2i - 1].x[2][2]
        x4[i] = sim4.sol.u[4i - 3].x[2][2]
    end
    if sim1.params.sf
        label_low = "sf: low - mid"
        label_mid = "sf: mid - high"
    else
        label_low = "low - mid"
        label_mid = "mid - high"
    end
    if invert
        xlabel!(p, "x")
        ylabel!(p, "t")
        plot!(p, x1 .- x2, sim1.sol.t; label=label_low)
        plot!(p, scale .* (x2 .- x4), sim1.sol.t; label=label_mid)
        # scatter!(p, x4, sim1.sol.t; label=label_high)
    else
        ylabel!(p, "x")
        xlabel!(p, "t")
        plot!(p, sim1.sol.t, x1 .- x2; label=label_low)
        plot!(p, sim1.sol.t, scale .* (x2 .- x4); label=label_mid)
    end
    savefig(p, "particle_location_errors.pdf")
    return p
end

function plot_position_conv_factor(sim1, sim2, sim4, scale=2; invert=false)
    pgfplotsx()
    p = plot(; title=L"Convergence Factor $= 2^p$")

    nt = length(sim1.sol.t)
    x1 = zeros(nt)
    x2 = zeros(nt)
    x4 = zeros(nt)
    for i in 1:nt
        x1[i] = sim1.sol.u[i].x[2][2]
        x2[i] = sim2.sol.u[2i - 1].x[2][2]
        x4[i] = sim4.sol.u[4i - 3].x[2][2]
    end
    err1 = x1 - x2
    err2 = x2 - x4
    co = err1 ./ err2

    xlabel!(p, "t")
    plot!(p, sim1.sol.t, co; label=false)
    savefig(p, "particle_location_convergence_factor.pdf")
    return p
end

end # end of module
