module Types

export SimulationParameters, SimulationParameters1D, SimulationParametersInteracting,
       Simulation,
       print_simulation_info

abstract type SimulationParameters{T} end

struct SimulationParameters1D{T} <: SimulationParameters{T}
    L::Int64
    N::Int64
    dx::T
    dt::T
    x::Vector{T}
    q::T
    m::T
    z0::T
    v0::T
    vmax::T
    A::T
    ω::T
    k::T
    direction::Int
    Nosc::T
    sf::Bool
    cfl::T
end

struct SimulationParametersInteracting{T} <: SimulationParameters{T}
    L::Int64
    N::Int64
    dx::T
    x::Vector{T}
    q1::T
    m1::T
    x1::T
    v1::T
    q2::T
    m2::T
    x2::T
    v2::T
    cfl::T
end

struct Simulation{T,S}
    params::SimulationParameters{T}
    sol::S
end

function print_simulation_info(title; params...)
    println("\n" * "="^50)
    println("Running Simulation: $title")
    println("-"^50)
    for (key, value) in params
        println(rpad(string(key), 20) * ": " * string(value))
    end
    println("="^50 * "\n")
end

function print_simulation_info(sim::SimulationParameters1D; title="")
    println("\n" * "="^50)
    if isempty(title)
        if sim.k == 0.0
            println("Running Simulation: Forced Motion")
            SimType = "Forced Motion"
        else
            println("Running Simulation: Particle in Potential")
            SimType = "Harmonic Potential"
        end
    else
        println("Running Simulation: $title")
        SimType = title
    end
    println("-"^50)

    T = 2π / sim.ω
    tf = sim.Nosc * T

    println(rpad("SimType", 20) * ": " * SimType)
    println(rpad("L", 20) * ": " * string(sim.L))
    println(rpad("N", 20) * ": " * string(sim.N))
    println(rpad("dx", 20) * ": " * string(sim.dx))
    println(rpad("dt", 20) * ": " * string(sim.dt))
    println(rpad("cfl", 20) * ": " * string(sim.cfl))
    println(rpad("Nosc", 20) * ": " * string(sim.Nosc))
    println(rpad("tf", 20) * ": " * string(tf))
    println(rpad("T", 20) * ": " * string(T))
    println(rpad("vmax", 20) * ": " * string(sim.vmax))
    println(rpad("q", 20) * ": " * string(sim.q))
    println(rpad("A", 20) * ": " * string(sim.A))
    println(rpad("omega", 20) * ": " * string(sim.ω))
    if sim.k != 0.0
        println(rpad("k", 20) * ": " * string(sim.k))
        println(rpad("m", 20) * ": " * string(sim.m))
    end

    println("="^50 * "\n")
end

end
