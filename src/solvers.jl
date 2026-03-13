



"""

    StepSolver <: AbstractStepSolver

Solve using a pre-defined number of steps


### Fields

+ `NT::Int` -- Number of steps
+ `savers::Vector{<:AbstractSavers}` -- List of AbstractSavers
+ `detector::AbstractDetector` -- Which detector to use
+ `samplers::Vector{<:AbstractSampler}` -- List of AbstractSamplers


### Notes

- The object may have different combinations of savers, detectors and samplers, not all must be set. If one is not set then it will be defined as empty and not used
- When detectors are used with the savers, then all AbstractSavers will be saved when the event ends.


### Examples

- `StepSolver(100)` -- Solve without detection, saving or sampling
- `StepSolver(100, detector::AbstractDetector)` -- Solve and detect events using an AbstractDetector
- `StepSolver(100, detector, savers=[saver1, saver2])` -- Solve, detect and save.
- `StepSolver(100, detector, samplers=[sampler1, sampler2])` -- Solve, detect and sample.
...


"""
struct StepSolver{S<:Vector{<:AbstractSaver}, D<:AbstractDetector, Sa<:Vector{<:AbstractSampler}} <: AbstractStepSolver
    NT::Int
    savers::S
    detector::D
    samplers::Sa

    function StepSolver(NT::Int, savers::S) where S
        detector = EmptyDetector()
        samplers = [EmptySampler()]
        new{typeof(savers), typeof(detector), typeof(samplers)}(NT, savers, detector, samplers)
    end

    function StepSolver(NT::Int, savers::S, detector::D) where {S, D<:AbstractDetector}
        samplers = [EmptySampler()]

        new{typeof(savers), typeof(detector), typeof(samplers)}(NT, savers, detector, samplers)
    end

    function StepSolver(NT::Int, savers::S, samplers::Sa) where {S, Sa<:Vector{AbstractSampler}}
        detector = EmptyDetector()

        new{typeof(savers), typeof(detector), typeof(samplers)}(NT, savers, detector, samplers)
    end

    function StepSolver(NT::Int, savers::S, detector::D, samplers::Sa) where {S, D, Sa}
        new{typeof(savers), typeof(detector), typeof(samplers)}(NT, savers, detector, samplers)
    end

    function StepSolver{S, D, Sa}(NT::Int, savers::S, detector::D, samplers::Sa) where {S, D, Sa}
        new{S, D, Sa}(NT, savers, detector, samplers)
    end
end

"""

    TimeSolver <: AbstractTimeSolver

Solve up to the specified time (in seconds).


### Fields

+ `tf::Int` -- Number of steps
+ `savers::Vector{<:AbstractSavers}` -- List of AbstractSavers
+ `detector::AbstractDetector` -- Which detector to use
+ `samplers::Vector{<:AbstractSampler}` -- List of AbstractSamplers


### Notes

- The object may have different combinations of savers, detectors and samplers, not all must be set. If one is not set then it will be defined as empty and not used
- When detectors are used with the savers, then all AbstractSavers will be saved when the event ends.


### Examples

- `TimeSolver(100)` -- Solve without detection, saving or sampling
- `TimeSolver(100, detector::AbstractDetector)` -- Solve and detect events using an AbstractDetector
- `TimeSolver(100, detector, savers=[saver1, saver2])` -- Solve, detect and save.
- `TimeSolver(100, detector, samplers=[sampler1, sampler2])` -- Solve, detect and sample.
...


"""
struct TimeSolver{S<:Vector{<:AbstractSaver}, D<:AbstractDetector, Sa<:Vector{<:AbstractSampler}} <: AbstractTimeSolver
    tf::Float64
    savers::S
    detector::D
    samplers::Sa

    function TimeSolver(tf::Float64, savers::S) where S
        detector = EmptyDetector()
        samplers = [EmptySampler()]
        new{typeof(savers), typeof(detector), typeof(samplers)}(tf, savers, detector, samplers)
    end

    function TimeSolver(tf::Float64, savers::S, detector::D) where {S, D<:AbstractDetector}
        samplers = [EmptySampler()]

        new{typeof(savers), typeof(detector), typeof(samplers)}(tf, savers, detector, samplers)
    end

    function TimeSolver(tf::Float64, savers::S, samplers::Sa) where {S, Sa<:Vector{<:AbstractSampler}}
        detector = EmptyDetector()

        new{typeof(savers), typeof(detector), typeof(samplers)}(tf, savers, detector, samplers)
    end

    function TimeSolver(tf::Float64, savers::S, detector::D, samplers::Sa) where {S, D, Sa}
        new{typeof(savers), typeof(detector), typeof(samplers)}(tf, savers, detector, samplers)
    end

    function TimeSolver{S, D, Sa}(tf::Float64, savers::S, detector::D, samplers::Sa) where {S, D, Sa}
        new{S, D, Sa}(tf, savers, detector, samplers)
    end
end


"""
    solve(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, solver::AbstractStepSolver)
    solve(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, solver::AbstractStepSolver, plotter::LivePlotter)
    solve(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, solver::AbstractTimeSolver)
    solve(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, solver::AbstractTimeSolver, plotter::LivePlotter)


Solve the problem in a limited number of steps or time with optional live plotting

### Arguments

- `experiment::AbstractExperiment` -- Experimental setup to solve
- `algorithm::AbstractAlgorithm` -- Which algorithm to use
- `solver::AbstractSolver` -- Which solver to use
- `plotter::LivePlotter` -- Live plot

### Notes

- When using GPUs live plotting is supported only with set_GPUbackend("unified")


"""
function solve(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, solver::AbstractStepSolver)

    state = experiment.state
    dx = state.dx
    V  = state.V
    theta  = state.theta

    NT = solver.NT

    detector = solver.detector
    savers = solver.savers
    samplers = solver.samplers



    for i in 1:NT

        dx, V, theta = algorithm(dx, V, theta)

        sample(samplers, stepper, state)
        detect(detector)
        simsave(savers)


    end

end

function solve(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, solver::AbstractStepSolver, plotter::LivePlotter)

    state = experiment.state
    dx = state.dx
    V  = state.V
    theta  = state.theta

    NT = solver.NT

    detector = solver.detector
    savers = solver.savers
    samplers = solver.samplers


    for i in 1:NT

        dx, V, theta = algorithm(dx, V, theta)

        sample(samplers, stepper, state)
        detect(detector)
        simsave(savers)
        UpdatePlot(plotter)


    end

end

function solve(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, solver::AbstractTimeSolver)

    state = experiment.state
    stepper = algorithm.stepper
    dx = state.dx
    V  = state.V
    theta  = state.theta

    tf = solver.tf

    detector = solver.detector
    savers = solver.savers
    samplers = solver.samplers




    while stepper.time <= tf

        dx, V, theta = algorithm(dx, V, theta)

        sample(samplers, stepper, state, detector.eventN)
        if typeof(detector) != EmptyDetector
            detect(detector, savers)
        else
            simsave(savers, stepper.step)
        end



    end

end

function solve(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, solver::AbstractTimeSolver, plotter::LivePlotter)

    state = experiment.state
    stepper = algorithm.stepper
    dx = state.dx
    V  = state.V
    theta  = state.theta

    tf = solver.tf

    detector = solver.detector
    savers = solver.savers
    samplers = solver.samplers


    while stepper.time <= tf

        dx, V, theta = algorithm(dx, V, theta)

        sample(samplers, stepper, state)
        detect(detector)
        simsave(savers)
        UpdatePlot(plotter)


    end
end



function benchmarksolve(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, steps, print_steps=false)

    state = experiment.state
    dx = state.dx
    V  = state.V
    theta  = state.theta


    for i in 1:steps
        if print_steps
            println("$(now()) step: $i")
        end
        dx, V, theta = algorithm(dx, V, theta)

    end


end
