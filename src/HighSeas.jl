module HighSeas

# using BenchmarkTools
using Makie
using AbstractFFTs
using LinearAlgebra
using FFTW
using FastBroadcast
import GPUArrays.AbstractGPUArray
using JLD2
using TimeZones
using StyledStrings
using CSV
using Glob
using SFTPClient


include("base.jl")
include("backends.jl")
include("materials.jl")
include("geometries.jl")
include("utils.jl")

include("experiments.jl")
include("laws/laws.jl")
include("steppers.jl")
include("algorithms.jl")
include("detectors.jl")
include("samplers.jl")
include("plotters.jl")
include("savers.jl")
include("solvers.jl")


const supported_GPU_platforms  = ["CUDA", "METAL"]
const supported_platforms = ["CPU", supported_GPU_platforms...]
const global_settings = Dict{String, Any}("backend"=>CPUBackend())

Base.copy(x::T) where T = T([deepcopy(getfield(x, k)) for k ∈ fieldnames(T)]...)

export global_settings
export State, Catalog
export SimpleMaterial
export Grid, PowerGrid, RectangleFault, RectanglePatch, CustomPatch, RectangleNucleation, Domain
export ReadSheet, RandomState, get_backend, get_available_platforms, get_available_GPUplatforms, set_CPUbackend, set_GPUbackend

export BP4QDExp
export StressFFT, AgeingLaw, ExplicitRate, LinearizedRate, HybridRate, GoverningEquations, DoubleError
export AdaptiveStepper
export CustomNewtonSolver
export SimpleDetector, CatalogDetector
export PointSampler, SectionSampler, ContourSampler, SamplerSaver
export RSPlotter
export StepSaver, CatalogSaver, SnaptshotSaver
export TimeSolver, StepSolver



end # module HighSeas
