
"""
    AgeingLaw <: AbstractStateLaw

Create an AgeingLaw object to use the ageing law for calculating the state in the rate and state simulation

### Fields

+ Dc::Float64 -- Dc value
+ theta_mem::AbstractArray{Float64} -- This is the memory matrix in which the result is saved

### Notes

- When using GPUs, it is possible to decide where the masks reside using `gpu_id`

### Examples

- `AgeingLaw(experiment::AbstractExperiment; gpu_id::Int=0)` -- AgeingLaw rate law for the current experiment
"""
struct AgeingLaw{M<:AbstractArray{Float64}} <: AbstractStateLaw

    Dc::Float64
    theta_mem::M # This is the memory matrix in which the result is saved

    function AgeingLaw(experiment::AbstractExperiment)

        Dc = experiment.material.Dc
        theta_mem = experiment.state.theta

        new{typeof(theta_mem)}(Dc, theta_mem)
    end

    function AgeingLaw{M}(Dc, theta_mem::M) where M
        new{M}(Dc, theta_mem)
    end


end

function (ageingLaw::AgeingLaw)(theta, V, dt)
    Dc = ageingLaw.Dc
    theta_mem = ageingLaw.theta_mem

    @.. thread=true theta_mem = theta * exp(-V * dt / Dc) + Dc / V * (1. - exp(-V * dt / Dc))

    return theta_mem
end


"""
SlipLaw(experiment::AbstractExperiment)

Create an SlipLaw object to use the slip law for calculating the state in the rate and state simulation

### Fields

+ `Dc::Float64` -- Dc value
+ `theta_mem::AbstractArray{Float64}` -- This is the memory matrix in which the result is saved

### Notes

- When using GPUs, it is possible to decide where the masks reside using `gpu_id`

### Examples

- `SlipLaw(experiment::AbstractExperiment; gpu_id::Int=0)` -- SlipLaw law for the current experiment
"""
struct SlipLaw <: AbstractStateLaw

    Dc::Float64
    theta_mem:: AbstractArray{Float64}

    function SlipLaw(experiment::AbstractExperiment)

        Dc = experiment.material.Dc
        theta = experiment.state.theta

        new(Dc, theta)
    end

end

function (slipLaw::SlipLaw)(theta, V, dt)
    Dc = ageingLaw.Dc
    theta_mem = ageingLaw.theta_mem

    @.. thread=true theta_mem = Dc/V * (V * theta / Dc)^exp(-V * dt / Dc)

    return theta_mem
end
