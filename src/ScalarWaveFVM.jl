module ScalarWaveFVM

include("Interpolations.jl")
include("InitialData.jl")
include("Equations.jl")
include("Limiters.jl")
include("Reconstructions.jl")
include("NumericalFluxes.jl")
include("ParticleMotion.jl")
include("ScalarField.jl")
# include("FractionalStepMethods.jl")
include("ODE.jl")
include("Run.jl")

include("Convergence.jl")
include("PlottingUtils.jl")
include("ForcedPlots.jl")
include("InteractingPlots.jl")
include("ParticlePotentialPlots.jl")
include("PreparePlots.jl")

end
