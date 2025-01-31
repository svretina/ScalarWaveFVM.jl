module Convergence

function convergence(L, N1, tf=1.0, cfl=0.5)
    x1, sol1 = ScalarWaveFVM.run(L, N1, tf, cfl)
    x2, sol2 = ScalarWaveFVM.run(L, 2N1, tf, cfl)
    # x4, sol4 = ScalarWaveFVM.run(L, 4N1, tf, cfl)
    if all(x1 .== x2[1:2:end]) #&& all(x1 .== x4[1:4:end])
        return x1, sol1, sol2# , sol4
    else
        throw("grids do not match")
    end
end

function plot_errors(ti, scale=2)
    t = sol1.t[ti]
    p = scatter(x, sol1[:, 1, ti] .- dtGaussian1D.(t, x, A, σ))
    scatter!(p, x, (sol2[1:2:end, 1, 2ti - 1] .- dtGaussian1D.(t, x, A, σ)) .* scale)
    return p
end

function plot_self_errors(ti, scale=2)
    plot1 = scatter(x, (sol2[1:2:end, 1, 2ti - 1] .- sol4[1:4:end, 1, 4ti - 3]) .* scale)
    scatter!(plot1, x, sol1[:, 1, ti] .- sol2[1:2:end, 1, 2ti - 1])
end

function plot_res(t)
    p = scatter(x, sol1[:, 1, t]; label="low res")
    scatter!(p, x, sol2[1:2:end, 1, 2t - 1]; label="mid res")
    #scatter!(p, x, sol4[1:4:end, 1, 4t - 3]; label="high res")
    xaxis!(p, "x")
    title!(p, "Π")
    return p
end

function plot_res_errors(t)
    p = plot_res(t)
    scatter!(p, x, sol1[:, 1, t] .- sol2[1:2:end, 1, 2t - 1]; label="low-mid")
    scatter!(p, x, (sol2[1:2:end, 1, 2t - 1] .- sol4[1:4:end, 1, 4t - 3]); label="mid-high")
end

function total_variation(q::AbstractArray)
    TV = zero(eltype(q))
    for i in eachindex(q)
        TV += abs(q[i] - q[i - 1])
    end
    return TV
end

function convergence_order(p)
    order = zeros(length(sol1.t))
    for i in 1:(length(sol1.t))
        t = sol1.t[i]
        tmp1 = sum(abs.(sol1[:, 1, i] .- dtGaussian1D.(t, x, A, σ)) .^ p)^(1 / p)
        tmp2 = sum(abs.(sol2[1:2:end, 1, 2i - 1] .- dtGaussian1D.(t, x, A, σ)) .^ p)^(1 / p)
        order[i] = log2(tmp1 / tmp2)
    end
    return order
end

function self_convergence_order(p)
    order = zeros(length(sol1.t))
    for i in 1:(length(sol1.t))
        tmp1 = sum(abs.(sol1[:, 1, i] .- sol2[1:2:end, 1, 2i - 1]) .^ p)^(1 / p)
        tmp2 = sum(abs.(sol2[1:2:end, 1, 2i - 1] .- sol4[1:4:end, 1, 4i - 3]) .^ p)^(1 /
                                                                                     p)
        order[i] = log2(tmp1 / tmp2)
    end
    return order
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

end # end of module
