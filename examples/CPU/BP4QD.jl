using HighSeas
using JLD2
using BenchmarkTools
# using CairoMakie
# using GLMakie

input_dict = readSheet("BP4input.txt")
# saved_state, step, time = load("BP4QD_out/2026-02-12T11:16:16.270/CPU/saved_StepSaver.jld2")["data"]

# Define the domain

grid = PowerGrid(input_dict)
fault = RectangleFault(input_dict, grid)

patch = RectanglePatch(input_dict, grid)
nucleation = RectangleNucleation(input_dict, grid)

domain = Domain(grid, fault, patch, nucleation)

# Define the material

material = SimpleMaterial(input_dict)


# Define the experiment
experiment = BP4QDExp(input_dict, material, domain, 10, "BP4QD_out")

# experiment = BP4QDExp(input_dict, material, domain, 10, saved_state)

# Define the Govenring equations

stresslaw = StressFFT(experiment)
statelaw = AgeingLaw(experiment)

explicitlaw = ExplicitRate(experiment)
linearlaw = LinearizedRate(experiment)
hybridlaw = HybridRate(0.00001, explicitlaw, linearlaw)

governing_equations = GoverningEquations(hybridlaw, stresslaw, statelaw)


# Define the error law and stepper

errorlaw = DoubleError(experiment)
stepper  = AdaptiveStepper(input_dict, errorlaw)
# stepper  = AdaptiveStepper(input_dict, errorlaw, saved_state)


# stepper = saved_step


# Define the algorithm used to solve

algorithm = CustomNewtonSolver(experiment, governing_equations, stepper)

# Define the detector (i.e. what to do when V reaches a certain value)

detector = CatalogDetector(1e-3, 1e-2, experiment, algorithm)

# Additional live plotting. This works only with CPU or GPU in UnifiedMemory, does not work with GPU in DeviceMemory

# plotter = RSPlotter(experiment, algorithm, 10)

# Additional samplers

samplers = Array{HighSeas.AbstractSampler, 1}(undef, 16)
for sp in 1:14
    samplers[sp] = PointSampler("BP4_sample_points.txt", sp, 700000, experiment)
end
samplers[15] = SectionSampler(0.0, "y", 700000, experiment)
samplers[16] = ContourSampler(:V, 1e-3, experiment)

# # Define savers

stepsaver = StepSaver(500, experiment, algorithm) # save simulation state at every 500 steps
catalogsaver = CatalogSaver(500, experiment, algorithm) # save catalog at every 500 steps

# snapshotsaver = SnaptshotSaver("BP4QD_out", 100, plotter, experiment, algorithm) # save a figure at every 500 steps
# CairoMakie.activate!()
samplersaver = SamplerSaver(500, samplers, experiment, algorithm) # save the samplers at every 500 steps

# savers = [stepsaver, catalogsaver]
# savers = [stepsaver, catalogsaver, snapshotsaver, samplersaver]
# or
# savers = [stepsaver, catalogsaver, snapshotsaver]
#or
savers = [stepsaver, catalogsaver, samplersaver]
#or ...




# Define the solver (i.e. use definite steps or time?)
tf = input_dict["tf"]*(365*24*60*60)

# solver = TimeSolver(tf, savers, detector)
# or add samplers
solver = TimeSolver(tf, savers, detector)

# Solve
HighSeas.solve(experiment, algorithm, solver)
# HighSeas.solve(experiment, algorithm, solver, plotter)
