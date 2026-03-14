using HighSeas
using JLD2
using CUDA
using Fractalizer
using GLMakie

set_GPUbackend() # To use UnifiedMemory add "unified" as the method argument. Default is DeviceMemory. AMDGPU ony supports DeviceMemory

input_dict = ReadSheet("BP4input.txt")
points = [[3e4, 3e4, -3e4, -3e4, 3e4] [1.5e4, -1.5e4, -1.5e4, 1.5e4, 1.5e4]] #RW patch points
np1 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[6442887322735277629, -8987213128142308954, -333252884332351366, 8464945807962482870])
np2 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[8513830690257299402, 5301462472252722888, 3588050925270478459, 3787793271851686014])
np3 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[6278251954719919174, 7730772275841400182, -2768929581059409545, 9035333999970238256])
np4 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[3947142509679233647, 626970696438225036, 7033574499115343085, -1445958421458504479])

templates = [random_template(n) for n in [np1,np2,np3,np4]]

shape = ClosedShape(points)
s = [1+2*input_dict["h"]/shape.l, 1+2*input_dict["h"]/shape.w]

buffer = shape * s

fractal = fractalize(shape, templates, input_dict["cellsizex"])

# Define the domain

grid = PowerGrid(input_dict)
fault = RectangleFault(input_dict, grid)

# patch = RectanglePatch(input_dict, grid)
patch = CustomPatch(fractal.points, grid)
input_dict["h"] = patch.h
nucleation = RectangleNucleation(input_dict, grid)

domain = Domain(grid, fault, patch, nucleation)


# Define the material

material = SimpleMaterial(input_dict)


# Define the experiment

experiment = BP4QDExp(input_dict, material, domain, 1000)

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
# or
savers = [stepsaver, catalogsaver, samplersaver]
#or ...




# Define the solver (i.e. use definite steps or time?)
tf = input_dict["tf"]*(365*24*60*60)

# solver = TimeSolver(tf, savers, detector)
# or add samplers
solver = TimeSolver(tf, savers, detector, samplers)

# Solve
HighSeas.solve(experiment, algorithm, solver)
