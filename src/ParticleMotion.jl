module ParticleMotion

function oscillator(t, x0, A, ω)
    arg = ω * t
    x = x0 + A * sin(arg)
    Aω = A * ω
    v = Aω * cos(arg)
    @assert abs(v) < 1.0
    a = -Aω * ω * sin(arg)
    return x, v, a
end

function oscillator1(t, x0, vmax, f0, direction)
    @assert vmax < 1.0
    ω = 2π * f0
    A = vmax / ω
    arg = ω * t
    prosimo = sign(direction)
    x = x0 + A * sin(prosimo * arg)
    v = prosimo * vmax * cos(arg)
    @assert abs(v) < 1.0
    a = -ω * prosimo * vmax * sin(arg)
    return x, v, a
end

function oscillator2(t, x0, vmax, A, direction)
    @assert vmax < 1.0
    f0 = vmax / (2π * A)
    ω = 2π * f0
    arg = ω * t
    prosimo = sign(direction)
    x = x0 + A * sin(prosimo * arg)
    v = prosimo * vmax * cos(arg)
    @assert abs(v) < 1.0
    a = -ω * prosimo * vmax * sin(arg)
    return x, v, a
end

end #end of module
