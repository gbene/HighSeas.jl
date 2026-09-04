include("rate.jl")
include("state.jl")
include("stress.jl")
include("errors.jl")

"""
    AbstractHybridLaw

Abstract type for hybrid laws that combine different laws based on some criterion.

"""
abstract type AbstractHybridLaw end

"""
    AbstractHybridRateLaw <: AbstractHybridLaw

Abstract type for hybrid rate laws that use different rate laws depending on a slip rate threshold.

"""
abstract type AbstractHybridRateLaw <: AbstractHybridLaw end

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