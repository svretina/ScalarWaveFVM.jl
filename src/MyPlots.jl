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

function make_gif_with_particle(x, sol, params)
    nt = length(sol.t)

    x10 = params.x1
    x20 = params.x2
    L = params.L
    direction1 = params.direction1
    direction2 = params.direction2

    anim = @animate for i in 1:10:nt
        p = plot()
        ylims!(p, -0.1, 2.5)
        t = sol.t[i]
        x1, v1, a1 = ParticleMotion.oscillator(t, x10, L, direction1)
        x2, v2, a2 = ParticleMotion.oscillator(t, x20, L, direction2)
        scatter!(p, x, sol[:, 1, i]; label=false)
        vline!(p, [x1]; label=false)
        vline!(p, [x2]; label=false)
        title!(p, "Π at t[$(i)]=$(round(t,digits=3))")
    end
    gif(anim, "anim_fps15.gif"; fps=15)
end

end # end of module
