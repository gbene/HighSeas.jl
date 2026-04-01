"""

Abstract type used for all algorithms used to solve an experiment

"""
abstract type AbstractAlgorithm end

"""

Abstract subtype used for Newton algorithm

"""
abstract type AbstractNewton <: AbstractAlgorithm end

"""

Abstract type used for all backends

"""
abstract type AbstractBackend end

"""

Abstract type used for GPU backends

"""
abstract type AbstractGPUBackend<:AbstractBackend end

"""

Abstract type used to define a simulation State object

"""
abstract type AbstractState end

"""

Abstract type used to define a simulation Catalog object

"""
abstract type AbstractCatalog end

"""

Abstract type used to define a Detector used to detect events

"""
abstract type AbstractDetector end

"""

Abstract type used to define an experiment

"""
abstract type AbstractExperiment end

"""

Abstract subtype used to define an benchmark experiment

"""
abstract type AbstractBenchExperiment <: AbstractExperiment end

"""

Abstract type used to define the underlying grids of the simulations

"""
abstract type AbstractGrid end

"""

Abstract subtypetype used to define a grid with power 2 elements and square cells

"""
abstract type AbstractPowerGrid <:AbstractGrid end

"""

Abstract type used to define a fault element in the simulation

"""
abstract type AbstractFault end

"""

Abstract type used to define a rate weakening element in the simulation

"""
abstract type AbstractPatch end


"""

Abstract type used to define a nucleation patch element in the simulation

"""
abstract type AbstractNucleation end


"""

Abstract type used to group all elements of the simulated geometry (i.e. fault, rate weakening patch etcetc).

### Notes

For now only one fault, patch and nucleation can be defined per domain. Also only one domain per experiment can be used.

"""
abstract type AbstractDomain end


"""

Abstract type used to define loaded data

"""
abstract type AbstractLoadedObject end


"""

Abstract type used to define a material used in the simulation

"""
abstract type AbstractMaterial end

"""

Abstract type used to define a plotter

"""
abstract type AbstractPlotter end


"""

Abstract subtype used to define a plotter that updates live while simulating

"""
abstract type LivePlotter <: AbstractPlotter end


"""

Abstract type used to define a sampler to gather data while the simulation runs

"""
abstract type AbstractSampler end


"""

Abstract type used to define savers

"""
abstract type AbstractSaver end


"""

Abstract type used to define a solver, i.e. how the AstractAlgorithm is managed

"""
abstract type AbstractSolver end


"""

Abstract subtype used to define solvers that run in a predefined number of steps

"""
abstract type AbstractStepSolver <: AbstractSolver end

"""

Abstract subtype used to define solvers that run in a predefined time

"""
abstract type AbstractTimeSolver <: AbstractSolver end

"""

Abstract type used to define the stepper of the algorithm, i.e. how dt is managed

"""
abstract type AbstractStepper end

"""

Abstract subtype used to define a stepper where dt adapts

"""
abstract type AbstractAdaptiveStepper <: AbstractStepper end


"""

Abstract type used to define laws that govern the simulation

"""
abstract type AbstractLaw end

"""

Abstract type used to group the rate-and-state laws

"""
abstract type AbstractGoverningEquations end

"""

Abstract subtype used to define State laws (i.e. aging or slip law)

"""
abstract type AbstractStateLaw <: AbstractLaw end
"""

Abstract subtype used to define Rate laws (i.e. linearized or explicit)

"""
abstract type AbstractRateLaw <: AbstractLaw end

"""

Abstract subtype used to define a hybrid law (use different laws depending on a threshold)

"""
abstract type AbstractHybridRateLaw <: AbstractRateLaw end

"""

Abstract subtype used to define Stress laws

"""
abstract type AbstractStressLaw <: AbstractLaw end

"""

Abstract subtype used to define Error laws used in AbstractAdaptiveSteppers

"""
abstract type AbstractErrorLaw <: AbstractLaw end
