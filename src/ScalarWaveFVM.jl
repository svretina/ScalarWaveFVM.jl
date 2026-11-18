module ScalarWaveFVM

export enable_plotting!, plotting_loaded

include("Interpolations.jl")
include("InitialData.jl")
include("Equations.jl")
include("Limiters.jl")
include("Reconstructions.jl")
include("NumericalFluxes.jl")
include("ParticleMotion.jl")
include("ScalarField.jl")
include("FractionalStepMethods.jl")
include("RetardedPotentials.jl")

include("ODE.jl")
include("Types.jl")
include("Run.jl")
include("Energies.jl")

const _PLOTTING_FILES = (
    "PlottingUtils.jl",
    "Convergence.jl",
    "ForcedPlots.jl",
    "InteractingPlots.jl",
    "ParticlePotentialPlots.jl",
    "EnergyBalanceIngredientPlots.jl",
    "PreparePlots.jl",
    "RetardedPotentialPlots.jl",
)

const _PLOTTING_MODULES = (
    :PlottingUtils,
    :Convergence,
    :ForcedPlots,
    :InteractingPlots,
    :ParticlePotentialPlots,
    :EnergyBalanceIngredientPlots,
    :PreparePlots,
    :RetardedPotentialPlots,
)

const _REVISE_PKGID = Base.PkgId(
    Base.UUID("295af30f-e4ad-537b-8983-00126c2a3abe"),
    "Revise",
)

"""
    plotting_loaded() -> Bool

Return `true` once the optional plotting modules have been loaded into
`ScalarWaveFVM`.
"""
plotting_loaded() = all(name -> isdefined(@__MODULE__, name), _PLOTTING_MODULES)

function _maybe_track_with_revise(path::AbstractString)
    revise = get(Base.loaded_modules, _REVISE_PKGID, nothing)
    revise === nothing && return nothing
    revise.track(@__MODULE__, path)
    return nothing
end

"""
    enable_plotting!()

Load the optional plotting submodules on demand. This keeps the core solver
loadable even when the active Julia environment cannot currently precompile the
Makie plotting stack.
"""
function enable_plotting!()
    plotting_loaded() && return nothing

    for file in _PLOTTING_FILES
        path = joinpath(@__DIR__, file)
        include(path)
        _maybe_track_with_revise(path)
    end

    return nothing
end

end
