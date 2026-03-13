abstract type AbstractLaw end
abstract type AbstractGoverningEquations end
abstract type AbstractStateLaw <: AbstractLaw end
abstract type AbstractRateLaw <: AbstractLaw end
abstract type AbstractHybridRateLaw <: AbstractRateLaw end
abstract type AbstractStressLaw <: AbstractLaw end
abstract type AbstractErrorLaw <: AbstractLaw end


include("rate.jl")
include("state.jl")
include("stress.jl")
include("errors.jl")

#Don't know if this struct is necessary or not.
# It can be handy and makes things more organized but at the same time it is an additional layer.

"""
    GoverningEquations <: AbstractGoverningEquation

Commodity object used to group the governing equations of the simulation

### Fields

+ ratelaw::AbstractRateLaw -- Rate law
+ stresslaw::AbstractStressLaw -- Stress law
+ statelaw::AbstractStateLaw -- State law



"""
struct GoverningEquations{RL<:AbstractRateLaw, TL<:AbstractStressLaw, SL<:AbstractStateLaw} <: AbstractGoverningEquations

    ratelaw::RL
    stresslaw::TL
    statelaw::SL

    function GoverningEquations(ratelaw::AbstractRateLaw, stresslaw::AbstractStressLaw, statelaw::AbstractStateLaw)
        new{typeof(ratelaw), typeof(stresslaw), typeof(statelaw)}(ratelaw, stresslaw, statelaw)
    end

end
