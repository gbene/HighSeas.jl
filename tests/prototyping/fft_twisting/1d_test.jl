using GLMakie
using FFTW
using Random


abstract type AbstractSlipFunc end

function phaseList(Q)
    theta = collect(2π .* (0:Q-1) ./ Q)
    theta[theta .> π] .-= 2π
    return theta
end

function fft_modes(N::Integer)
    if iseven(N)
        return collect(vcat(0:(div(N,2)-1), -div(N,2):-1))
    else
        return collect(vcat(0:div(N-1,2), -div(N-1,2):-1))
    end
end

struct CircleSlipFunc <: AbstractSlipFunc
    Lx::Float64
    r::Float64
end

function (slip::CircleSlipFunc)(x)
    
    s = slip.r^2 - x^2

    # display(mask)

    if s > 0
        return sqrt(s)
    else
        return 0.0
    end
end

struct PhaseShift
    N::Float64
    θ::Float64
    
end


function (phase_shift::PhaseShift)(x)
    # return exp(-1im * phase_shift.θ)
    return exp(-1im * phase_shift.θ*x/phase_shift.N)
    # return exp(-1im * phase_shift.θ/phase_shift.L)
end

function plotphaseAvg(x, theta_angles, delta, Q, Lx)

    tau = zeros(length(delta))
    fig = Figure(font=20)
    ax = Axis(fig[2,1], aspect = DataAspect(), title="1D Slip", ylabel="δ(x)")
    axend = Axis(fig[2, 5], title="τ, phaseavg", xlabel="x", ylabel="τ")

    lines!(ax, x, delta)

    Nx = length(delta)
    kx0_vec = (2π / Lx) .* fft_modes(Nx)
    ix = collect(0:Nx-1)

    G = 3464^2*2670 #cs^2*rho
    nu = 0.25
    for (i, theta) in enumerate(theta_angles)
        
        axtw = Axis(fig[i,2], title="1D Slip, twisted θ = $(rad2deg(theta))", xlabel="x", ylabel="δ(x)")
        # axA = Axis(fig[i,3], title="1D Slip, twisted θ = $(rad2deg(theta))", xlabel="x", ylabel="δ(x)")

        axuntw = Axis(fig[i,4], title="1D Slip, untwisted θ = $(rad2deg(theta))", xlabel="x", ylabel="τ")
        display(theta)
        kx_vec = kx0_vec .+ theta / Lx


        display(kx_vec)
        A = -G/2 * (kx_vec/(1-nu))
        # A = kx_vec

        # lines!(axA, x, A)

        # shifter = PhaseShift(Nx, theta) 
        # tw = shifter.(ix)  
        tw = exp.(-im .* theta .* ix ./ Nx)
        untw = conj.(tw)
        delta_shift = tw .* delta
        
        lines!(axtw, x, real(delta_shift))
        lines!(axtw, x, delta, alpha=0.6)

        delta_hat = fft(delta_shift)

        tau_hat = delta_hat #.* kx_vec

        tau_theta = ifft(tau_hat).* untw

        lines!(axuntw, x, real(tau_theta))

        tau += tau_theta

    end

    tau /= Q
    lines!(axend, x, real(tau))
    display(fig)
    # return tau
end


Lx = 10.0

r = 5
x = -Lx/2:0.01:Lx/2

Nx = length(x)

slip_func = CircleSlipFunc(Lx, r)
shift0_func = PhaseShift(Nx, 0.0)
shiftpi_func = PhaseShift(Nx, π)

twist_theta0 = shift0_func.(x)
untwist_theta0 = conj(twist_theta0)

twist_thetapi = shiftpi_func.(x)
untwist_thetapi = conj(twist_thetapi)


delta = slip_func.(x)
Q = 2
theta_angles = phaseList(Q)


plotphaseAvg(x, theta_angles, delta, Q, Lx)

# fig = Figure(fontsize = 20)



# delta_shift0 = twist_theta0 .* delta
# delta_shiftpi = twist_thetapi .* delta


# delta_hat0 = fft(delta_shift0)
# delta_hatpi = fft(delta_shiftpi)

# tau_hat0 = delta_hat0 .* rand(Xoshiro(0))
# tau_hatpi = delta_hatpi .* rand(Xoshiro(4))

# tau_0 = ifft(tau_hat0).* untwist_theta0
# tau_pi = ifft(tau_hatpi).* untwist_thetapi



# tau = (tau_0+tau_pi)/2


# ax = Axis(fig[2,1], aspect = DataAspect(), title="1D Slip", ylabel="δ(x)")


# ax2 = Axis(fig[1,2], aspect = DataAspect(), title="1D Slip, twisted θ=0", ylabel="δ(x)")
# ax3 = Axis(fig[2,2], aspect = DataAspect(), title="1D slip, twisted θ=π", xlabel="x", ylabel="δ(x)")

# # ax4 = Axis(fig[1,3], title="ℱ[δ(x)], twisted θ=0", xlabel="x", ylabel="ℱ[δ(x)] (real)")
# # ax5 = Axis(fig[3,3], title="ℱ[δ(x)], twisted θ=π", xlabel="x", ylabel="ℱ[δ(x)] (real)")

# ax6 = Axis(fig[1,4], aspect = DataAspect(), title="τ, untwisted θ=0", xlabel="x", ylabel="τ")
# ax7 = Axis(fig[2,4], aspect = DataAspect(), title="τ, untwisted θ=π", xlabel="x", ylabel="τ")
# ax8 = Axis(fig[2,5], aspect = DataAspect(), title="τ, phaseavg", xlabel="x", ylabel="τ")




# lines!(ax, x, delta)
# ylims!(ax, -r*1.1, r*1.1)

# lines!(ax2, x, real(delta_shift0))
# # lines!(ax2, x, delta, alpha=0.6)
# ylims!(ax2, -r*1.1, r*1.1)

# lines!(ax3, x, real(delta_shiftpi))
# # lines!(ax3, x, delta, alpha=0.6)
# ylims!(ax3, -r*1.1, r*1.1)

# # lines!(ax4, x, real.(fftshift(delta_hat0)))
# # lines!(ax5, x, real.(fftshift(delta_hatpi)))


# lines!(ax6, x, real.(tau_0))
# lines!(ax7, x, real.(tau_pi))
# ylims!(ax6, -r*1.1, r*1.1)
# ylims!(ax7, -r*1.1, r*1.1)


# lines!(ax8, x, real.(tau))
# ylims!(ax8, -r*1.1, r*1.1)



# # rowsize!(fig.layout, 1, Aspect(1, 1))

# # lines!(ax5, x, real.(delta_hatpi))


# display(fig)