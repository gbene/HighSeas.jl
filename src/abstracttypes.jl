abstract type AbstractAlgorithm end
abstract type AbstractNewton <: AbstractAlgorithm end

abstract type AbstractBackend end
abstract type AbstractGPUBackend<:AbstractBackend end

abstract type AbstractState end
abstract type AbstractCatalog end

abstract type AbstractDetector end

abstract type AbstractExperiment end
abstract type AbstractBenchExperiment <: AbstractExperiment end

abstract type AbstractGrid end
abstract type AbstractPowerGrid <:AbstractGrid end

abstract type AbstractFault end
abstract type AbstractPatch end
abstract type AbstractNucleation end

abstract type AbstractDomain end

abstract type AbstractLoadedObject end

abstract type AbstractMaterial end

abstract type AbstractPlotter end
abstract type LivePlotter <: AbstractPlotter end

abstract type AbstractSampler end

abstract type AbstractSaver end

abstract type AbstractSolver end
abstract type AbstractStepSolver <: AbstractSolver end
abstract type AbstractTimeSolver <: AbstractSolver end

abstract type AbstractStepper end
abstract type AbstractAdaptiveStepper <: AbstractStepper end
