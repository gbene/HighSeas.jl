abstract type AbstractMaterial end


"""
SimpleMaterial(input_dict::Dict)

Create a SimpleMaterial

# Attributes


+ cs::Float64       # Shear wave speed
+ rho::Float64      # Density
+ nu::Float64       # Poissons ratio
+ a::Float64        # a value of the rate weakening part
+ aRs::Float64      # a value of the rate strengthening part
+ b::Float64        # b value of rate and state
+ fr::Float64       # Reference friction coefficient
+ G::Float64        # Youngs modulus
+ eta::Float64      # Radiation damping coefficient
+ Dc::Float64       # Characteristic state evolution distance

# Notes
"""
struct SimpleMaterial <: AbstractMaterial
    cs::Float64
    rho::Float64
    nu::Float64
    a::Float64
    aRs::Float64
    b::Float64
    fr::Float64
    G::Float64
    eta::Float64
    Dc::Float64

    function SimpleMaterial(input_dict::Dict)
        cs  = input_dict["cs"]
        rho = input_dict["rho"]
        nu  = input_dict["nu"]
        a   = input_dict["a"]
        aRS = input_dict["aRS"]
        b   = input_dict["b"]
        fr  = input_dict["fr"]
        Dc  = input_dict["Dc"]
        G   = cs^2*rho
        eta = G/(2*cs)

        new(cs, rho, nu, a, aRS, b, fr, G, eta, Dc)
    end

    function SimpleMaterial(cs, rho, nu, a, aRS, b, fr, G, eta, Dc)
        new(cs, rho, nu, a, aRS, b, fr, G, eta, Dc)
    end
end
