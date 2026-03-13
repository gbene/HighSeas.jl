

struct EmptySampler <: AbstractSampler
end

"""
    PointSampler <: AbstractSampler

Sample simulation at a given point

### Fields

- mask::AbstractArray{Int8} -- Preallocated array used to define the point in the domain (1 only at the coordinate)
- NT::Int -- Number of samples
- dxs::Vector{Float64} -- Sampled values of slip (m)
- Vs::Vector{Float64} -- Sampled values of slip rate (m/s)
- thetas::Vector{Float64} -- Sampled values of θ
- taus::Vector{Float64} -- Sampled values of τ (Pa)
- times::Vector{Float64} -- Sampled times
- sample_point_id::Int -- Id of the sample point
- sample_point_x::Float64 -- x coord of the sample point (m)
- sample_point_y::Float64 -- y coord of the sample point (m)

### Notes

- When using GPUs, it is possible to decide where the mask reside using `gpu_id`
- It is possible to use a sample point file structured as follows

```
n x y
1 -36.0 0.0
2 -22.5 -7.50
3 -16.5 -12.0
...
```


### Examples

-- `PointSampler(sample_point_id::Int, sample_point_x::Float64, sample_point_y::Float64, NT::Int, experiment::AbstractExperiment; gpu_id::Int=0)` -- Sample at point
-- `PointSampler(sample_points_paths::String, sample_point_id, NT, experiment; gpu_id=0)` -- Sample point at a poisition given by a txt file

"""
struct PointSampler{B<:AbstractArray{Int8}} <: AbstractSampler

    # state::S
    # stepper::ST
    # temp::M
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


    function PointSampler(sample_points_paths::String, sample_point_id::Int, NT::Int, experiment::AbstractExperiment; gpu_id::Int=0)

        # state = experiment.state
        # stepper = algorithm.stepper

        sample_points = CSV.File(open(sample_points_paths))

        grid = experiment.domain.grid
        sample_point_x = sample_points.x[sample_point_id]*1000
        sample_point_y = sample_points.y[sample_point_id]*1000

        mask = @. Int8((grid.x == sample_point_x) * (grid.y == sample_point_y));
        # temp = zeros(grid.n_elementsy, grid.n_elementsx)

        dxs           = fill(NaN, NT)
        Vs            = fill(NaN, NT)
        thetas        = fill(NaN, NT)
        taus          = fill(NaN, NT)
        times         = fill(NaN, NT)

        if typeof(get_backend()) <: AbstractGPUBackend
            mask = memcopy(mask, gpu_id)
            # temp = memcopy(temp)
        end

        new{typeof(mask)}(mask, NT, dxs, Vs, thetas, taus, times, sample_point_id, sample_point_x, sample_point_y)
    end

    function PointSampler(sample_point_id::Int, sample_point_x::Float64, sample_point_y::Float64, NT::Int, experiment::AbstractExperiment; gpu_id::Int=0)

        grid = experiment.domain.grid

        mask = @. Int8((grid.x == sample_point_x) * (grid.y == sample_point_y));
        # temp = zeros(grid.n_elementsy, grid.n_elementsx)

        dxs           = fill(NaN, NT)
        Vs            = fill(NaN, NT)
        thetas        = fill(NaN, NT)
        taus          = fill(NaN, NT)
        times         = fill(NaN, NT)

        if typeof(get_backend()) <: AbstractGPUBackend
            mask = memcopy(mask, gpu_id)
        end

        new{typeof(mask)}(mask, NT, dxs, Vs, thetas, taus, times, sample_point_id, sample_point_x, sample_point_y)
    end

    function PointSampler{B}(mask::B, NT, dxs, Vs, thetas, taus, times, sample_point_id, sample_point_x, sample_point_y) where {B}
        new{B}(mask, NT, dxs, Vs, thetas, taus, times, sample_point_id, sample_point_x, sample_point_y)
    end

end


function (pointSampler::PointSampler)(stepper, state)

    mask = pointSampler.mask
    # temp = pointSampler.temp
    step = stepper.step

    dx = state.dx
    V = state.V
    theta = state.theta
    tau = state.tau
    t = stepper.time



    dxs = pointSampler.dxs
    Vs = pointSampler.Vs
    thetas = pointSampler.thetas
    taus   = pointSampler.taus
    times = pointSampler.times

    # @.. thread=true temp = dx * mask

    dxs[step]       = dot(dx, mask)

    # @.. thread=true temp = V * mask

    Vs[step]        = dot(V, mask)

    # @.. thread=true temp = theta * mask

    thetas[step]    = dot(theta, mask)

    # @.. thread=true temp = tau * mask

    taus[step]      = dot(tau, mask)

    times[step]     = t

    return nothing

end



"""
    SectionSampler <: AbstractSampler

Sample simulation using a section

### Fields

- `grid::AbstractGrid` -- Grid useed in the simulation
- `temp::AbstractArray{Float64}` -- Temporary matrix used to calculate the section values
- `mask::AbstractArray{Int8}` -- Temporary mask used to calculate the section values
- `NT::Int` -- Number of section to sample
- `quantity::Symbol` -- Which quantity of AbstractState to sample
- `section::Matrix{Float64}` -- NT x (n+1) matrix where n is the grid resolution and NT is the number of section to sample
- `coord::Float64` -- Position to sample
- `axis::String` -- Which axis

### Notes

- When using GPUs, it is possible to decide where the mask reside using `gpu_id`
- When using GPUs `temp` will always be set in unified memory
- The first column of the matrix defines the maximum slip (useful to define seismic/aseismic slips)


### Examples

-- `SectionSampler("dx", coord=0.0, axis="y", NT::Int, experiment::AbstractExperiment; gpu_id::Int=0)` -- Sample slip at all points where y=0.0
"""
struct SectionSampler{G<:AbstractGrid, M<:AbstractArray{Float64}, B<:AbstractArray{Int8}} <: AbstractSampler

    grid::G
    temp::M
    mask::B

    NT::Int

    quantity::Symbol
    section::Matrix{Float64}
    coord::Float64
    axis::String


    function SectionSampler(quantity::String, coord::Float64, axis::String, NT::Int, experiment::AbstractExperiment; gpu_id::Int=0)

        # state = experiment.state
        # stepper = algorithm.stepper


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
            mask = memcopy(mask, gpu_id)
            old_memtype = HighSeas.get_backend().memtype # I really don't like this but it is necessary as section is on CPU memory.
            HighSeas.set_GPUbackend("unified")

            temp = memcopy(temp, gpu_id)

            HighSeas.set_GPUbackend(old_memtype)

        end

        new{typeof(grid), typeof(temp), typeof(mask)}(grid, temp, mask, NT, Symbol(quantity), section, coord, axis)
    end

    function SectionSampler{G, M, B}(grid::G, temp::M, mask::B, NT, quantity, section, coord, axis) where {G,M,B}
        new{G, M, B}(grid, temp, mask, NT, quantity, section, coord, axis)
    end


end

function (sectionSampler::SectionSampler)(stepper, state)

    mask = sectionSampler.mask
    temp = sectionSampler.temp
    step = stepper.step

    V = state.V
    data = getproperty(state, sectionSampler.quantity)


    section = sectionSampler.section

    @.. thread=true temp = data * mask


    section[step,1]       = maximum(V)
    if sectionSampler.axis == "y"
        section[step,2:end]   = sum(temp, dims=1)
    else
        section[step,2:end]   = sum(temp, dims=2)
    end

    return nothing

end


"""
    ContourSampler <: AbstractSampler

Sample rupture contours

### Fields

+ `contour::AbstractArray{Float64}` -- Contour matrix
+ `mask::AbstractArray{Int8}` -- Temporary mask used to calculate the rupture contours
+ `thresh::Float64` -- Threshold value for the contour
+ `first_contour::Int8` -- Pre allocated matrix used to store the first contour
+ `t_fc::Float64` -- First contour time
+ `quantity::Symbol` -- Quantity to sample

### Notes

- When using GPUs, it is possible to decide where the mask reside using `gpu_id`


### Examples

-- `ContourSampler("V", 1e-3, experiment::AbstractExperiment; gpu_id=0)` -- Sample slip speed when reaching 1e-3
"""
mutable struct ContourSampler{M<:AbstractArray{Float64}, B<:AbstractArray{Int8}} <: AbstractSampler

    contour::M
    mask::B


    thresh::Float64
    first_contour::Int8
    t_fc::Float64
    quantity::Symbol

    function ContourSampler(quantity::String, thresh::Float64, experiment::AbstractExperiment; gpu_id::Int=0)

        state = experiment.state
        # stepper = algorithm.stepper


        sz = size(experiment.domain.grid.x)

        # temp = zeros(grid.n_elementsy, grid.n_elementsx)
        contour = fill(NaN, sz)

        mask = zeros(Int8, sz)


        first_contour = 0
        t_fc = 0.0

        if typeof(get_backend()) <: AbstractGPUBackend
            contour = memcopy(contour, gpu_id)
            mask    = memcopy(mask, gpu_id)
            # temp    = memcopy(temp)
        end

        quantity = Sumbol(quantity)

        new{typeof(contour), typeof(mask)}(contour, mask, thresh, first_contour, t_fc, quantity)
    end

    function ContourSampler{M, B}(contour::M, mask::B, thresh, first_contour, t_fc, quantity) where {M, B}
        new{M, B}(contour, mask, thresh, first_contour, t_fc, quantity)
    end


end

function (contourSampler::ContourSampler)(stepper, state, eventN)


    # eventN = contourSampler.detector.eventN


    if eventN == 1
        mask = contourSampler.mask
        t = stepper.time

        field = getfield(state, contourSampler.quantity)
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



function sample(sampler::AbstractSampler, stepper, state)

    if typeof(sampler) != EmptySampler
        sampler(stepper, state)
    end

    return nothing
end

function sample(sampler::ContourSampler, stepper, state, eventN::Int)

    if typeof(sampler) != EmptySampler
        sampler(stepper, state, eventN)
    end

    return nothing
end

function sample(samplers::Vector{<:AbstractSampler}, stepper, state, eventN)
    for sampler in samplers
        sample(sampler, stepper, state, eventN)
    end
    return nothing
end
