abstract type AbstractSampler end


struct EmptySampler <: AbstractSampler
end


struct PointSampler{S<:AbstractState, ST<:AbstractStepper, M<:AbstractArray{Float64}, B<:AbstractArray{Int8}} <: AbstractSampler

    state::S
    stepper::ST
    temp::M
    mask::B

    NT::Int


    dxs::Vector{Float64}
    Vs::Vector{Float64}
    thetas::Vector{Float64}
    taus::Vector{Float64}
    times::Vector{Float64}

    sample_point_id::Int
    sample_point_x::Float64
    sample_point_y::Float64


    function PointSampler(sample_points_paths::String, sample_point_id::Int, NT::Int, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

        state = experiment.state
        stepper = algorithm.stepper

        sample_points = CSV.File(open(sample_points_paths))

        grid = experiment.domain.grid
        sample_point_x = sample_points.x[sample_point_id]*1000
        sample_point_y = sample_points.y[sample_point_id]*1000

        mask = @. Int8((grid.x == sample_point_x) * (grid.y == sample_point_y));
        temp = zeros(grid.n_elementsy, grid.n_elementsx)

        dxs           = fill(NaN, NT)
        Vs            = fill(NaN, NT)
        thetas        = fill(NaN, NT)
        taus          = fill(NaN, NT)
        times         = fill(NaN, NT)

        if typeof(get_backend()) <: AbstractGPUBackend
            mask = memcopy(mask)
            temp = memcopy(temp)
        end

        new{typeof(state), typeof(stepper), typeof(temp), typeof(mask)}(state, stepper, temp, mask, NT, dxs, Vs, thetas, taus, times, sample_point_id, sample_point_x, sample_point_y)
    end

    function PointSampler{S,ST,M,B}(state::S, stepper::ST, temp::M, mask::B, NT, dxs, Vs, thetas, taus, times, sample_point_id, sample_point_x, sample_point_y) where {S,ST,M,B}
        new{S,ST,M,B}(state, stepper, temp, mask, NT, dxs, Vs, thetas, taus, times, sample_point_id, sample_point_x, sample_point_y)
    end

end


function (pointSampler::PointSampler)()

    mask = pointSampler.mask
    temp = pointSampler.temp
    step = pointSampler.stepper.step

    dx = pointSampler.state.dx
    V = pointSampler.state.V
    theta = pointSampler.state.theta
    tau = pointSampler.state.tau
    t = pointSampler.stepper.time



    dxs = pointSampler.dxs
    Vs = pointSampler.Vs
    thetas = pointSampler.thetas
    taus   = pointSampler.taus
    times = pointSampler.times

    @.. thread=true temp = dx * mask

    # display(sum(temp))

    dxs[step]       = sum(temp)

    @.. thread=true temp = V * mask

    Vs[step]        = sum(temp)


    @.. thread=true temp = theta * mask

    thetas[step]        = sum(temp)

    @.. thread=true temp = tau * mask

    taus[step]      = sum(temp)

    times[step]     = t

    return nothing

end


struct SectionSampler{S<:AbstractState, ST<:AbstractStepper,G<:AbstractGrid, M<:AbstractArray{Float64}, B<:AbstractArray{Int8}} <: AbstractSampler

    state::S
    stepper::ST
    grid::G
    temp::M
    mask::B

    NT::Int


    section::Matrix{Float64}
    coord::Float64
    axis::String



    function SectionSampler(coord::Float64, axis::String, NT::Int, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

        state = experiment.state
        stepper = algorithm.stepper


        grid = experiment.domain.grid

        if axis == "x"
            mask = @. Int8(grid.x == coord*1000)
        elseif axis == "y"
            mask = @. Int8(grid.y == coord*1000)
        else
            error("direction not supported")
        end
        temp = zeros(grid.n_elementsy, grid.n_elementsx)

        section = fill(NaN, (NT, sum(mask)+1))

        if typeof(get_backend()) <: AbstractGPUBackend
            mask = memcopy(mask)
            old_memtype = HighSeas.get_backend().memtype
            HighSeas.set_GPUbackend("unified")

            temp = memcopy(temp)

            HighSeas.set_GPUbackend(old_memtype)

        end

        new{typeof(state), typeof(stepper), typeof(grid), typeof(temp), typeof(mask)}(state, stepper, grid, temp, mask, NT, section, coord, axis)
    end


end

function (sectionSampler::SectionSampler)()

    mask = sectionSampler.mask
    temp = sectionSampler.temp
    step = sectionSampler.stepper.step

    V = sectionSampler.state.V
    dx = sectionSampler.state.dx


    section = sectionSampler.section

    @.. thread=true temp = dx * mask


    section[step,1]       = maximum(V)
    section[step,2:end]   = sum(temp, dims=1)

    return nothing

end

mutable struct ContourSampler{S<:AbstractState, ST<:AbstractStepper, D<:AbstractDetector, M<:AbstractArray{Float64}, B<:AbstractArray{Int8}} <: AbstractSampler

    state::S
    stepper::ST
    detector::D
    temp::M
    mask::B


    contour::M
    thresh::Float64
    first_contour::Int8
    t_fc::Float64
    field::M

    function ContourSampler(field::Symbol, thresh::Float64, detector::AbstractDetector, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

        state = experiment.state
        stepper = algorithm.stepper


        grid = experiment.domain.grid

        temp = zeros(grid.n_elementsy, grid.n_elementsx)

        mask = zeros(Int8, size(temp))

        contour = fill(NaN, size(temp))

        first_contour = 0
        t_fc = 0.0

        if typeof(get_backend()) <: AbstractGPUBackend
            mask    = memcopy(mask)
            temp    = memcopy(temp)
            contour = memcopy(contour)
        end

        field = getproperty(state, field)

        new{typeof(state), typeof(stepper), typeof(detector), typeof(temp), typeof(mask)}(state, stepper, detector, temp, mask, contour, thresh, first_contour, t_fc, field)
    end


end

function (contourSampler::ContourSampler)()


    eventN = contourSampler.detector.eventN



    if eventN == 1
        mask = contourSampler.mask
        temp = contourSampler.temp
        step = contourSampler.stepper.step
        t = contourSampler.stepper.time
        state = contourSampler.state

        field = contourSampler.field
        contour = contourSampler.contour


        maxF = maximum(field)

        @.. thread=true mask = (field >= contourSampler.thresh) * (isnan(contour))
        if (maxF >= contourSampler.thresh) && contourSampler.first_contour == 0
            contourSampler.first_contour = 1;
            contourSampler.t_fc = t;
        end
        @.. thread=true contour = iszero(mask) * contour + mask * (t-contourSampler.t_fc);

    end

    return nothing

end



function sample(sampler::AbstractSampler)

    if typeof(sampler) != EmptySampler
        sampler()
    end

    return nothing
end

function sample(samplers::Vector{<:AbstractSampler})
    for sampler in samplers
        sample(sampler)
    end
    return nothing
end
