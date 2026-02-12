abstract type AbstractSolver end
abstract type AbstractStepSolver <: AbstractSolver end
abstract type AbstractTimeSolver <: AbstractSolver end




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

        sample(samplers)
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

        sample(samplers)
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

        sample(samplers)
        detect(detector)
        simsave(savers)



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

        sample(samplers)
        detect(detector)
        simsave(savers)
        UpdatePlot(plotter)


    end
end



function benchmarksolve(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, steps)

    state = experiment.state
    dx = state.dx
    V  = state.V
    theta  = state.theta


    for i in 1:steps

        dx, V, theta = algorithm(dx, V, theta)

    end


end
