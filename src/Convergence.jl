module Convergence

function convergence(L, N1, tf=1.0, cfl=0.5, sf=true, muscl=true)
    x1, sol1 = ScalarWaveFVM.Run.run(L, N1, tf, cfl, sf, muscl)
    x2, sol2 = ScalarWaveFVM.Run.run(L, 2N1, tf, cfl, sf, muscl)
    x4, sol4 = ScalarWaveFVM.Run.run(L, 4N1, tf, cfl, sf, muscl)
    return x1, x2, x4, sol1, sol2, sol4
end

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
