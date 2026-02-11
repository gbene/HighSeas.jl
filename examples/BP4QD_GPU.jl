using HighSeas
using CUDA # Change this to AMDGPU, METAL and so on to use those backends. METAL does not support Float64 so for now it is not supported


set_GPUbackend() # To use UnifiedMemory add "unified" as the method argument. Default is DeviceMemory. AMDGPU ony supports DeviceMemory
input_dict = ReadSheet("BP4input.txt")

# Define the domain

grid = PowerGrid(input_dict)
fault = RectangleFault(input_dict, grid)

patch = RectanglePatch(input_dict, grid)
nucleation = RectangleNucleation(input_dict, grid)

domain = Domain(grid, fault, patch, nucleation)

# Define the material

material = SimpleMaterial(input_dict)


# Define the experiment

experiment = BP4QDExp(input_dict, material, domain, 12)

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



# Define the algorithm used to solve

algorithm = CustomNewtonSolver(experiment, governing_equations, stepper)

# Define the detector (i.e. what to do when V reaches a certain value)

detector = CatalogDetector(1e-3, 1e-2, experiment, algorithm)

# Additional live plotting. This works only with CUDA.UnifiedMemory

# plotter = RSPlotter(experiment, algorithm, 10)

# Additional samplers

# samplers = Array{HighSeas.AbstractSampler, 1}(undef, 15)
# for sp in 1:14
#     samplers[sp] = PointSampler("BP4_sample_points.txt", sp, 70000, experiment, algorithm)
# end
# samplers[15] = SectionSampler(0.0, "y", 70000, experiment, algorithm)

# Define savers

stepsaver = StepSaver("BP4QD_out", 500, experiment, algorithm) # save simulation state at every 500 steps
catalogsaver = CatalogSaver("BP4QD_out", 500, experiment, algorithm) # save catalog at every 500 steps

# snapshotsaver = SnaptshotSaver("BP4QD_out", 100, plotter, experiment, algorithm) # save a figure at every 500 steps
# samplersaver = SamplerSaver("BP4QD_out", 500, samplers, experiment, algorithm) # save the samplers at every 500 steps


savers = [stepsaver, catalogsaver]
# or
# savers = [stepsaver, catalogsaver, snapshotsaver]
#or
# savers = [stepsaver, catalogsaver, samplersaver]
#or ...




# Define the solver (i.e. use definite steps or time?)
tf = input_dict["tf"]*(365*24*60*60)

solver = TimeSolver(tf, savers, detector, samplers)

# Solve
HighSeas.solve(experiment, algorithm, solver)
