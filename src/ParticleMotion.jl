module ParticleMotion

# We fix the motion to be an oscillating
# particle in the domain [-L,L]
# starting at x0=x0 with a0=0 and v0=0.5
function oscillator(t, x0, L, direction)
    #f0 is the frequency that gives v=1
    _L = one(L) / L
    # f0 = 1/ (2π * L/4)
    # f0 = 2.0 * _L / π
    # ω = 2π * f0 / 2
    ω = 2.0 * _L
    arg = ω * t
    prosimo = sign(direction)
    x = x0 + 0.25L * sin(prosimo*arg)
    v = prosimo*0.5cos(arg)  
    a = -prosimo*sin(arg) * _L
    return (x, v, a)
end

end #end of module
