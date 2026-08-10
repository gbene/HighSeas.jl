using HighSeas
using CUDA # Change this to AMDGPU, METAL and so on to use those backends. METAL does not support Float64 so for now it is not supported
using JLD2
using GLMakie
using BenchmarkTools


set_GPUbackend() # To use UnifiedMemory add "unified" as the method argument. Default is DeviceMemory. AMDGPU ony supports DeviceMemory
input_dict = readSheet("circular_input.txt")
# saved_step = loadData("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/Asperity_out/2026_05_15T13_03_01.128/CUDA/saved_StepSaver.jld2")
# Define the domain

grid = PowerGrid(input_dict)
fault = HighSeas.CircularFault(input_dict, grid, 0.8)
patch = HighSeas.FractalAsperity(input_dict, grid, 0.9, 0.35,16)

domain = Domain(grid, fault, patch)

material = SimpleMaterial(input_dict)


# Define the experiment

experiment = HighSeas.Experiment(input_dict, material, domain, 100000, "Asperity_out")
# experiment = HighSeas.Experiment(input_dict, material, domain, 100000, "Asperity_out", saved_step)



stresslaw = StressFFT(experiment)
statelaw = AgeingLaw(experiment)

explicitlaw = ExplicitRate(experiment)
linearlaw = LinearizedRate(experiment)
hybridlaw = HybridRate(0.00001, explicitlaw, linearlaw)

governing_equations = GoverningEquations(hybridlaw, stresslaw, statelaw)


# Define the error law and stepper

errorlaw = DoubleError(experiment)
stepper  = AdaptiveStepper(input_dict, errorlaw)
# stepper  = AdaptiveStepper(input_dict, errorlaw, saved_step)<



# Define the algorithm used to solve

algorithm = CustomNewtonSolver(experiment, governing_equations, stepper)

# Define the detector (i.e. what to do when V reaches a certain value)

detector = CatalogDetector(1e-3, 1e-2, experiment, algorithm)

# Additional live plotting. This works only with CUDA.UnifiedMemory

# plotter = RSPlotter(experiment, algorithm, 30)

# Define savers

stepsaver = StepSaver(500, experiment, algorithm) # save simulation state at every 500 steps
catalogsaver = CatalogSaver(500, experiment, algorithm) # save catalog at every 500 steps


savers = [stepsaver, catalogsaver]



# Define the solver (i.e. use definite steps or time?)
tf = input_dict["tf"]*(365*24*60*60)

# solver = StepSolver(1, savers, detector)
solver = TimeSolver(tf, savers, detector)

# Solve
HighSeas.solve(experiment, algorithm, solver)
