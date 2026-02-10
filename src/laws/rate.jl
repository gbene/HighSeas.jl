
"""
ExplicitRate(experiment)

Create an ExplicitRate object to calculate explicitly the rate from the Rate and State equation.

# Attributes

eta::Float64                    # Radiation damping coefficient
fr::Float64                     # Friction reference coefficient
Vr::Float64                     # Reference slip velocity
Dc::Float64                     # Characteristic length
V_mem::AbstractArray{Float64}   # This is the memory matrix in which the result is saved


# Notes
"""
struct ExplicitRate{M<:AbstractArray{Float64}} <: AbstractRateLaw

    eta::Float64
    fr::Float64
    Vr::Float64
    Dc::Float64
    V_mem::M

    function ExplicitRate(experiment::AbstractExperiment)
        material = experiment.material
        eta = material.eta
        fr  = material.fr
        Dc  = material.Dc
        Vr  = experiment.Vr

        V_mem = experiment.state.V
        new{typeof(V_mem)}(eta, fr, Vr, Dc, V_mem)
    end

    function ExplicitRate{M}(eta, fr, Vr, Dc, V_mem) where M
        new{M}(eta, fr, Vr, Dc, V_mem)
    end

end

"""
LinearizedRate(experiment)

Create an LinearizedRate object to calculate the rate using a first order Taylor
expansion.

# Attributes

eta::Float64                    # Radiation damping coefficient
fr::Float64                     # Friction reference coefficient
Vr::Float64                     # Reference slip velocity
Dc::Float64                     # Characteristic length
V_mem::AbstractArray{Float64}   # This is the memory matrix in which the result is saved


# Notes
"""
struct LinearizedRate{M<:AbstractArray{Float64}} <: AbstractRateLaw

    eta::Float64
    fr::Float64
    Vr::Float64
    Dc::Float64
    V_mem::M

    function LinearizedRate(experiment::AbstractExperiment)
        material = experiment.material
        eta = material.eta
        fr  = material.fr
        Dc  = material.Dc
        Vr  = experiment.Vr

        V_mem = experiment.state.V

        new{typeof(V_mem)}(eta, fr, Vr, Dc, V_mem)
    end

    function LinearizedRate{M}(eta, fr, Vr, Dc, V_mem) where M
        new{typeof(V_mem)}(eta, fr, Vr, Dc, V_mem)
    end
end

"""
HybridRate(experiment)

Create an HybridRate object to calculate the rate using either the explicit or linearized method.

# Attributes

+ threshold::Float64            # Velocity threshold to decide which law to use
+ explicit::AbstractRateLaw     # First law to use in case the V is <= than the thresh
+ linear::AbstractRateLaw       # Second law to use in case the V is > than the thresh


# Notes
"""
struct HybridRate{E<:AbstractRateLaw, L<:AbstractRateLaw} <: AbstractHybridRateLaw

    threshold::Float64
    explicit::E
    linear::L

    function HybridRate(threshold::Float64, explicit, linear)
        new{typeof(explicit), typeof(linear)}(threshold, explicit, linear)
    end

    function HybridRate{E, L}(threshold::Float64, explicit::E, linear::L) where {E, L}
        new{typeof(explicit), typeof(linear)}(threshold, explicit, linear)
    end

end


function (explicitRate::ExplicitRate)(V, tau, theta, si, a, b)

    eta = explicitRate.eta
    fr  = explicitRate.fr
    Vr  = explicitRate.Vr
    Dc   = explicitRate.Dc
    V_mem = explicitRate.V_mem

    # This is more legible than a single string

    function calculateRate(V, tau, theta, si, a, b)

        alpha = (tau - eta * V) / (si * a);
        beta  = (fr + b * log(Vr * theta / Dc)) / a;
        return Vr * (exp(alpha - beta) - exp(-alpha - beta));

    end

    @.. thread=true V_mem = calculateRate(V, tau, theta, si, a, b)

    return V_mem

end
function (linearizedRate::LinearizedRate)(V, tau, theta, si, a, b)

    eta = linearizedRate.eta
    fr  = linearizedRate.fr
    Vr  = linearizedRate.Vr
    Dc   = linearizedRate.Dc
    V_mem = linearizedRate.V_mem

    # This is more legible than a single string

    function calculateRate(V, tau, theta, si, a, b)

        alpha = (tau - eta * V) / (si * a);
        beta  = (fr + b * log(Vr * theta / Dc)) / a;
        BB    = eta / si / a;

        dV = (Vr * (exp(alpha - beta) - exp(-alpha - beta)) - V)/(BB * Vr * (exp(alpha - beta) + exp(-alpha - beta)) + 1.);

        return V+dV

    end

    @.. thread=true V_mem = calculateRate(V, tau, theta, si, a, b)

    return V_mem

end

function (hybridRate::HybridRate)(V, tau, theta, si, a, b)
    maxV = maximum(V)
    explicit = hybridRate.explicit
    linear = hybridRate.linear

    if maxV > hybridRate.threshold
        V_mem = linear(V, tau, theta, si, a, b)
    else
        V_mem = explicit(V, tau, theta, si, a, b)
    end
    return V_mem
end
