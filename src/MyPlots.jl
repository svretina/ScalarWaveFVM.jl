module MyPlots

using Plots
using ..ParticleMotion

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

function plot_res(i)
    t1 = sol1.t[i]
    t2 = sol2.t[2i - 1]
    t4 = sol4.t[4i - 3]
    @assert t1 == t2
    @assert t2 == t4
    t = t1
    p = scatter(x1, sol1[:, 1, i]; label="low res")
    scatter!(p, x1, mymean(sol2[:, 1, 2i - 1]); label="mid res")
    # plot!(p, x1[begin]:0.01:x1[end], true_sol.(t, x1[begin]:0.01:x1[end]))
    scatter!(p, x1, mymean(mymean(sol4[:, 1, 4i - 3])); label="high res")
    xaxis!(p, "x")
    yaxis!(p, "Π")
    title!(p, "t=$(t)")
    return p
end

function plot_res_errors(i, true_sol)
    p = plot_res(i, true_sol)
    t = sol1.t[i]
    scatter!(p, x1, sol1[:, 1, i] .- true_sol.(t, x1); label="low-mid")
    scatter!(p, x1, mean(sol2[:, 1, 2i - 1]) .- true_sol.(t, x1); label="mid-high")
    return p
end

function make_gif_errors(scale=1.0)
    nt = length(sol1.t)
    p = plot()
    ylims!(p, -0.1, 1.1)
    anim = @animate for i in 1:10:nt
        plot_errors(i, scale)
    end
    gif(anim, "scaled_errors_fps15.gif"; fps=15)
end

function make_gif(x, sol)
    nt = length(sol.t)
    p = plot()
    ylims!(p, -0.1, 1.1)
    anim = @animate for i in 1:10:nt
        scatter(p, x, sol[:, 1, i])
    end
    gif(anim, "anim_fps15.gif"; fps=15)
end

function make_gif_with_forced_particle(x, sol, params)
    nt = length(sol.t)

    x10 = params.x1
    x20 = params.x2
    L = params.L
    direction1 = params.direction1
    direction2 = params.direction2
    vmax = params.vmax
    f0 = params.f0
    A = params.A
    anim = @animate for i in 1:10:nt
        p = plot()
        ylims!(p, -3.1, 5.5)
        t = sol.t[i]
        x1, v1, a1 = ParticleMotion.oscillator2(t, x10, vmax, A, direction1)
        x2, v2, a2 = ParticleMotion.oscillator2(t, x20, vmax, A, direction2)
        plot!(p, x, sol[:, 1, i]; label=false)
        vline!(p, [x1]; label=false)
        vline!(p, [x2]; label=false)
        title!(p, "Π at t[$(i)]=$(round(t,digits=3))")
    end
    gif(anim, "anim_fps15.gif"; fps=15)
end

function make_gif_with_particle(xs, sol)
    nt = length(sol.t)
    ymins = zeros(nt)
    ymaxs = zeros(nt)
    for i in 1:nt
        ymins[i] = min(sol.u[i].x[1][:, 1]...)
        ymaxs[i] = max(sol.u[i].x[1][:, 1]...)
    end
    ymin = min(ymins...)
    ymax = max(ymaxs...)
    p = plot()
    anim = @animate for i in 1:10:nt
        p = plot()
        # ylims!(p, 1.1ymin, 1.1ymax)
        t = sol.t[i]
        x1 = sol.u[i].x[2][2]

        plot!(p, xs, sol.u[i].x[1][:, 1]; label=false)
        vline!(p, [x1]; label=false)
        title!(p, "Π at t[$(i)]=$(round(t,digits=3))")
    end
    gif(anim, "anim_fps15.gif"; fps=15)
end

function make_gif_with_particles(x, sol)
    nt = length(sol.t)
    ymins = zeros(nt)
    ymaxs = zeros(nt)
    for i in 1:nt
        ymins[i] = min(sol.u[i].x[1][:, 1]...)
        ymaxs[i] = max(sol.u[i].x[1][:, 1]...)
    end
    ymin = min(ymins...)
    ymax = max(ymaxs...)
    anim = @animate for i in 1:10:nt
        p = plot()
        # ylims!(p, 1.1ymin, 1.1ymax)
        t = sol.t[i]
        x1 = sol.u[i].x[2][2]
        x2 = sol.u[i].x[2][5]

        plot!(p, x, sol.u[i].x[1][:, 1]; label=false)
        vline!(p, [x1]; label=false)
        vline!(p, [x2]; label=false)
        title!(p, "Π at t[$(i)]=$(round(t,digits=3))")
    end
    gif(anim, "anim_fps15.gif"; fps=15)
end

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

function plot_position1(sol; invert=false)
    p = plot()

    nt = length(sol.t)
    x1 = zeros(nt)
    x2 = zeros(nt)
    for (i, ti) in enumerate(sol.t)
        x1[i] = sol.u[i].x[2][2]
    end
    if invert
        plot!(p, x1, sol.t; title="particle 1", xlabel="x", ylabel="t", label=false)
    else
        plot!(p, sol.t, x1; title="particle 1", xlabel="t", ylabel="x", label=false)
    end
    savefig(p, "particle_location.png")
    return p
end

end # end of module
