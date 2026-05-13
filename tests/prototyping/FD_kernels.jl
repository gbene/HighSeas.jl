using SpecialFunctions
using PCHIPInterpolation
using Prony
using GLMakie

"""
    cumtrapz(X, Y)
    cumtrapz(X, F)
    cumtrapz!(Z, X, Y)
    cumtrapz!(Z, X, F)



Calculate the integral of a 1D function using the trapezoidal rule

### Arguments

    - `X::AbstractVector` -- Input vector
    - `Y::AbstractVector` -- Input vector
    - `F::Function` -- Function used to integrate over X
    - `Z::AbstractVector` -- Output vector

### Output

    - `AbstractVector` of the integral at the given point i

### Notes

Taken from: https://stackoverflow.com/questions/58139195/cumulative-integration-options-with-julia
"""
function cumtrapz(X, Y::T) where {T <: AbstractVector}
  # Check matching vector length
  @assert length(X) == length(Y)
  # Initialize Output
  out = similar(X)
  out[1] = 0
  # Iterate over arrays
  for i in 2:length(X)
    out[i] = out[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
  end
  # Return output
  out
end
function cumtrapz(X::T, F::Function) where {T <: AbstractVector}
    out = similar(X)
    out[1] = 0
    for i in 2:length(X)
        x_val = X[i]
        x_val_past = X[i-1]
        out[i] = out[i-1] + 0.5*(x_val - x_val_past)*(F(x_val) + F(x_val_past))
  end
  # Return output
  out

end


function cumtrapz!(Z::T, X::T, Y::T) where {T <: AbstractVector}
  # Check matching vector length
  @assert length(X) == length(Y)
  @assert length(Z) == length(X)
  # Initialize Output
  Z[1] = 0
  # Iterate over arrays
  for i in 2:length(X)
    Z[i] = Z[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
  end
  # Return output
  out
end
function cumtrapz!(Z::T, X::T, F::Function) where {T <: AbstractVector}
  # Check matching vector length
  @assert length(Z) == length(X)
  # Initialize Output
  Z[1] = 0
  # Iterate over arrays
  for i in 2:length(X)
    x_val = X[i]
    x_val_past = X[i-1]
    Z[i] = Z[i-1] + 0.5*(x_val - x_val_past)*(F(x_val) + F(x_val_past))
  end
  # Return output
  out
end

"""
    C2(rho)

Calculate the CII term as A2 in Lapusta and Liu 2009

### Arguments

    - `rho::AbstractVector` -- vector of ρ values
    - `gamma::AbstractVector` -- Vector of cₚ/cₛ * ρ values

### Output

`AbstractVector` of C2 values

### Bibliography

Lapusta, Nadia, and Yi Liu. ‘Three-Dimensional Boundary Integral Modeling of Spontaneous Earthquake Sequences and Aseismic Slip’.
Journal of Geophysical Research: Solid Earth 114, no. B9 (2009). https://doi.org/10.1029/2008JB005934.
"""
function C2(rho::AbstractVector{T}, gamma:: AbstractVector{T}, cp_cs::T) where {T}
    C2 = C3.(rho) + 4 .* rho .* (W(gamma) - W(rho)) - 4 ./ cp_cs .* besselj0.(gamma) + 3 .* besselj0.(rho)

    return C2
end

function C2(rho:: AbstractVector{T}, gamma:: AbstractVector{T}, cp_cs::T, W_interp::Interpolator) where {T}
    C2 = C3.(rho) + 4 .* rho .* (W_interp.(gamma) - W_interp.(rho)) - 4 ./ cp_cs .* besselj0.(gamma) + 3 .* besselj0.(rho)

    return C2
end


"""
    C3(rho)

Calculate the CIII term as A3 in Lapusta and Liu 2009

### Bibliography

Lapusta, Nadia, and Yi Liu. ‘Three-Dimensional Boundary Integral Modeling of Spontaneous Earthquake Sequences and Aseismic Slip’.
Journal of Geophysical Research: Solid Earth 114, no. B9 (2009). https://doi.org/10.1029/2008JB005934.
"""
function C3(rho::T) where T
    if rho > 0.0
        return besselj1(rho)/rho
    else
        return 0.5 #lim_{x->0} f(x)/x = 1/2
    end

end


"""
    W(rho)

Calculate the W term given ρ as A4 in Lapusta and Liu 2009

### Bibliography

Lapusta, Nadia, and Yi Liu. ‘Three-Dimensional Boundary Integral Modeling of Spontaneous Earthquake Sequences and Aseismic Slip’.
Journal of Geophysical Research: Solid Earth 114, no. B9 (2009). https://doi.org/10.1029/2008JB005934.
"""
function W(rho::T) where {T <: AbstractVector}
    1 .- cumtrapz(rho, C3)
end


"""
    lapustakernels(nu, rhoMax, npoints)

Calculate the convolution kernels from Lapusta and Liu 2009

### Arguments

    - `nu::T` -- Poisson ratio value
    - `rhoMax::T` -- Maximum value of the dimensionless Kernel parameter ρ
    - `dx::T` -- Resolution of the kernel defined as step
    - `npoints::Int` -- Resolution of the kernel defined as number of points

### Output
    - `rho::StepRange` -- Rho values as a range
    - `K2::AbstractVector` -- Vector of KII values
    - `K3::AbstractVector` -- Vector of KIII values

### Bibliography

Lapusta, Nadia, and Yi Liu. ‘Three-Dimensional Boundary Integral Modeling of Spontaneous Earthquake Sequences and Aseismic Slip’.
Journal of Geophysical Research: Solid Earth 114, no. B9 (2009). https://doi.org/10.1029/2008JB005934.

"""
function lapustakernels(nu::T, rhoMax::T, dx::T) where {T}
    cp_cs = sqrt(2*(1-nu)/(1-2*nu))

    rhoW = 0:dx:cp_cs*rhoMax
    rho  = 0:dx:rhoMax
    gamma = cp_cs .* rho


    W_vals = W(rhoW)

    W_interp   = Interpolator(rhoW, W_vals)

    C2_interp = C2(rho, gamma, cp_cs, W_interp)
    C2_vals   = C2(rho, gamma, cp_cs)

    K2_interp = 2*(1-1/cp_cs^2) .- cumtrapz(rho, C2_interp)
    K2 = 2*(1-1/cp_cs^2) .- cumtrapz(rho, C2_vals)
    K3_interp = W_interp.(rho)
    K3 = W(rho)


    length(K2)
    return rho, K2, K2_interp, K3, K3_interp


end

rho, K2, K2int, K3, K3int = lapustakernels(0.25, 180.0, 1e-3)




KC_fig = Figure()

K2_ax = Axis(KC_fig[1,1], title="K2 vs K2 cubic", xlabel="ρ", ylabel="KII")

K3_ax = Axis(KC_fig[1,2], title="K3 vs K3 cubic", xlabel="ρ", ylabel="KIII")

K2diff_ax = Axis(KC_fig[2,1], title="Difference", xlabel="ρ", ylabel="KII-KIIc")
K3diff_ax = Axis(KC_fig[2,2], title="Difference", xlabel="ρ", ylabel="KIII-KIIIc")

L = Label(KC_fig[0,:], "Kernels: Formulation v.s. Cubic", fontsize=30)

lines!(K2_ax, rho, K2, label="KII", linewidth=3)
lines!(K2_ax, rho, K2int, label="KII, cubic", linestyle=:dash, linewidth=3)
lines!(K3_ax, rho, K3, label="KIII", linewidth=3)
lines!(K3_ax, rho, K3int, label="KIII, cubic", linestyle=:dash, linewidth=3)
lines!(K2diff_ax, rho, K2 .- K2int)
lines!(K3diff_ax, rho, K3 .- K3int)
axislegend(K2_ax)
axislegend(K3_ax)

xlims!(K2_ax, 0,20)
xlims!(K3_ax, 0,20)



KP_fig = Figure()

K2P_ax = Axis(KP_fig[1,1], title="K2 vs Prony", xlabel="ρ", ylabel="KII")
K3P_ax = Axis(KP_fig[1,2], title="K3 vs Prony", xlabel="ρ", ylabel="KIII")


L = Label(KP_fig[0,:], "Kernels: Formulation v.s. Prony", fontsize=30)
KC_fig
