# Native Julia SOE fit for the Lapusta-Liu KII and KIII kernels.
#
# This is a direct Julia translation of the MATLAB prototype
# `fit_Lapusta_kernels_SOE.m`: a grouped greedy dictionary fit using real
# decaying exponentials and damped sine/cosine conjugate-pair groups.
#
# Run from the package root, for example:
#     julia --project=. examples/resources/fit_lapusta_kernels_soe.jl
#
# Output:
#     examples/resources/kernel_fits/lapusta_kernel_fit.jld2
#
# The JLD2 file contains both the compact SOE pole/residue data and the dense
# exact kernel table needed by the optional brute-force validation method.

using LinearAlgebra
using Statistics
using SpecialFunctions
using JLD2

struct FitSettings
    nu::Float64
    rhoWindow::Float64
    zeroExtensionFactor::Float64
    rhoLinearEnd::Float64
    nLinear::Int
    nLog::Int
    nTail::Int
    drAux::Float64
    taperFraction::Float64
    decayReal::Vector{Float64}
    decayPair::Vector{Float64}
    omegaPair::Vector{Float64}
    nGroupsKII::Int
    nGroupsKIII::Int
    weightFloor::Float64
    weightPower::Float64
    extraZeroWeight::Float64
end

function default_settings()
    return FitSettings(
        0.25,      # nu
        200.0,     # rhoWindow
        1.8,       # zeroExtensionFactor, what is this?
        10.0,      # rhoLinearEnd,this is used for indicating a linear distr at the start (until rho 10)
        300,       # nLinear, number of points for the first linear part
        500,       # nLog, number of points for the exp/log part
        220,       # nTail, number of points for the last linear part
        1e-3,      # drAux what is this? Seems to be the step for rho
        0.20,      # taperFraction
        collect(exp.(range(log(5e-2), log(5.0), length=12))), # decay real range
        collect(exp.(range(log(5e-2), log(5.0), length=14))), # decay pair range
        collect(range(0.4, 1.6, length=21)), # Omega pairs range
        8,         # KII selected groups
        4,         # KIII selected groups
        0.02,      # weight foor
        0.25,      # weight power
        2.0,       # extraZeroWeight
    )
end

function cumtrapz_uniform(y, dx)
    out = similar(y)
    out[1] = zero(eltype(y))
    @inbounds for i in 2:length(y)
        out[i] = out[i-1] + 0.5 * dx * (y[i] + y[i-1])
    end
    return out
end

function lininterp(x::AbstractVector, y::AbstractVector, xq::AbstractVector)
    out = similar(xq, Float64)
    n = length(x)
    @inbounds for iq in eachindex(xq)
        z = xq[iq]
        if z <= x[1]
            out[iq] = y[1]
        elseif z >= x[end]
            out[iq] = y[end]
        else
            i = searchsortedlast(x, z)
            t = (z - x[i]) / (x[i+1] - x[i])
            out[iq] = (1-t)*y[i] + t*y[i+1]
        end
    end
    return out
end

function compute_exact_kernels(nu, rhoMax, drAux)
    cp_cs = sqrt(2*(1-nu)/(1-2*nu))

    rhoW = collect(0.0:drAux:(cp_cs*rhoMax))
    J1_over_r = similar(rhoW)
    J1_over_r[1] = 0.5
    @inbounds for i in 2:length(rhoW)
        J1_over_r[i] = besselj1(rhoW[i]) / rhoW[i]
    end
    F = cumtrapz_uniform(J1_over_r, drAux)
    W = 1 .- F

    rho = collect(0.0:drAux:rhoMax)
    W_rho = lininterp(rhoW, W, rho)
    W_gamma = lininterp(rhoW, W, cp_cs .* rho)

    J1r = similar(rho)
    J1r[1] = 0.5
    @inbounds for i in 2:length(rho)
        J1r[i] = besselj1(rho[i]) / rho[i]
    end

    CII = J1r .+ 4 .* rho .* (W_gamma .- W_rho) .- 4 ./ cp_cs .* besselj0.(cp_cs .* rho) .+ 3 .* besselj0.(rho)
    KII = 2 * (1 - 1 / cp_cs^2) .- cumtrapz_uniform(CII, drAux)
    KIII = W_rho

    return (rho=rho, KII=KII, KIII=KIII, cp_cs=cp_cs)
end

function build_fit_grid(settings::FitSettings)
    rhoWindow = settings.rhoWindow
    rhoExtended = settings.zeroExtensionFactor * rhoWindow
    rhoLinEnd = min(settings.rhoLinearEnd, rhoWindow)
    rhoLinear = collect(range(0.0, rhoLinEnd, length=settings.nLinear))

    rhoLog = Float64[]
    if rhoWindow > rhoLinEnd * (1 + 1e-12)
        rhoLogStart = max(rhoLinEnd + 1e-6, 1e-3)
        rhoLog = collect(exp.(range(log(rhoLogStart), log(rhoWindow), length=settings.nLog)))
    end

    rhoTail = Float64[]
    if rhoExtended > rhoWindow
        drTail = (rhoExtended - rhoWindow) / max(settings.nTail, 1)
        rhoTail = collect((rhoWindow + drTail):drTail:rhoExtended)
    end

    return sort(unique(vcat(rhoLinear, rhoLog, rhoTail)))
end

function tail_taper_window(rho, rhoWindow, taperFraction)
    w = zeros(length(rho))
    taperFraction <= 0 && return Float64.(rho .<= rhoWindow)
    taperStart = (1 - taperFraction) * rhoWindow
    for i in eachindex(rho)
        if rho[i] <= taperStart
            w[i] = 1.0
        elseif rho[i] <= rhoWindow
            s = (rho[i] - taperStart) / (rhoWindow - taperStart)
            w[i] = 0.5 * (1 + cos(π*s))
        else
            w[i] = 0.0
        end
    end
    return w
end

function build_group_dictionary(rho, settings::FitSettings)
    groupsmat = zeros(length(rho), length(settings.decayReal)+length(settings.decayPair)*length(settings.omegaPair))
    groups = Vector{Matrix{Float64}}()
    meta = Vector{NamedTuple}()

    for (i, a) in enumerate(settings.decayReal)
        groupsmat[:,i] = exp.(-a .* rho)
        push!(groups, reshape(exp.(-a .* rho), :, 1))
        push!(meta, (kind=:real, a=a, omega=0.0))
    end
    last_index = length(settings.decayReal)
    for a in settings.decayPair
        E = exp.(-a .* rho)
        for omega in settings.omegaPair
            # display(last_index)
            groupsmat[:,last_index]   = E .* cos.(omega .* rho)
            groupsmat[:,last_index+1] = E .* sin.(omega .* rho)

            push!(groups, hcat(E .* cos.(omega .* rho), E .* sin.(omega .* rho)))
            push!(meta, (kind=:pair, a=a, omega=omega))
            last_index+=1
        end
    end

    return groups, groupsmat, meta
end

function select_groups_greedy(groups, target, weights, nGroups)

    nCandidates = length(groups)
    available = trues(nCandidates)
    selected = Int[]
    yW = weights .* target
    Bsel = zeros(length(target), 0)
    coef = zeros(0)
    display(min(nGroups, nCandidates))
    for _ in 1:min(nGroups, nCandidates)
        bestErr = Inf
        bestIdx = 0
        bestB = Bsel
        bestCoef = coef
        for idx in findall(available)
            Bcand = hcat(Bsel, groups[idx])
            Bw = Bcand .* weights
            coefCand = Bw \ yW
            residual = target - Bcand * coefCand
            err = norm(weights .* residual)
            if err < bestErr
                bestErr = err
                bestIdx = idx
                bestB = Bcand
                bestCoef = coefCand
            end
        end
        push!(selected, bestIdx)
        available[bestIdx] = false
        Bsel = bestB
        coef = bestCoef
        @info "selected group" iteration=length(selected) group=bestIdx weighted_error=bestErr
    end

    Bw = Bsel .* weights
    coef = Bw \ yW
    return selected, coef, Bsel
end

function groups_to_poles_residues(selected, coef, meta)
    poles = ComplexF64[]
    residues = ComplexF64[]
    idx = 1
    for g in selected
        m = meta[g]
        if m.kind == :real
            push!(poles, complex(-m.a, 0.0))
            push!(residues, complex(coef[idx], 0.0))
            idx += 1
        else
            A = coef[idx]
            B = coef[idx+1]
            c = 0.5 * complex(A, -B)
            p = complex(-m.a, m.omega)
            push!(poles, p)
            push!(poles, conj(p))
            push!(residues, c)
            push!(residues, conj(c))
            idx += 2
        end
    end
    return poles, residues
end

function fit_kernel(name, rhoFit, target, groups, meta, nGroups, settings; tapered=true)

    scale = maximum(abs.(target))

    scale = scale <= 0 ? 1.0 : scale

    weights = 1.0 ./ (settings.weightFloor*scale .+ abs.(target)).^settings.weightPower

    if tapered
        weights[rhoFit .> settings.rhoWindow] .*= settings.extraZeroWeight # this makes no sense as when tapered the values after rhoWindow are 0
    end

    selected, coef, Bsel = select_groups_greedy(groups, target, weights, nGroups)
    display(selected)

    # poles, residues = groups_to_poles_residues(selected, coef, meta)

    # yfit = evaluate_soe(rhoFit, poles, residues)
    # relL2 = norm(target - yfit) / max(norm(target), eps())
    # relInf = maximum(abs.(target-yfit)) / max(maximum(abs.(target)), eps())
    # @info "fit complete" name=name relL2=relL2 relInf=relInf npoles=length(poles)
    return selected, coef, Bsel
end

function evaluate_soe(rho, poles, residues)
    y = zeros(ComplexF64, length(rho))
    for (p,c) in zip(poles, residues)
        @. y += c * exp(p * rho)
    end
    return real.(y)
end

settings = default_settings()
rhoExtended = settings.zeroExtensionFactor * settings.rhoWindow

exact = compute_exact_kernels(settings.nu, rhoExtended, settings.drAux)
rhoFit = build_fit_grid(settings)
taper = tail_taper_window(rhoFit, settings.rhoWindow, settings.taperFraction)


KIIbase = lininterp(exact.rho, exact.KII, rhoFit)
KIIIbase = lininterp(exact.rho, exact.KIII, rhoFit)
KIItarget = KIIbase .* taper
KIIItarget = KIIIbase .* taper
# KIItarget[rhoFit .> settings.rhoWindow] .= 0.0 # this is unnecessary as the taper is already 0 after rhoWindow
# KIIItarget[rhoFit .> settings.rhoWindow] .= 0.0 # this is unnecessary as the taper is already 0 after rhoWindow


groups, groupsmat, meta = build_group_dictionary(rhoFit, settings)

# @info "fitting KII"
selected, coef, Bsel = fit_kernel("KII tapered", rhoFit, KIItarget, groups, meta,
    settings.nGroupsKII, settings; tapered=true)


idxls, cls, els = lsomp(groupsmat, KIItarget; invert = true, verbose = true, ϵrel = eps(), maxterms=8)


# @info "fitting KIII"
# polesIII, residuesIII = fit_kernel("KIII tapered", rhoFit, KIIItarget, groups, meta,
#     settings.nGroupsKIII, settings; tapered=true)

# outdir = joinpath(@__DIR__, "kernel_fits")
# mkpath(outdir)
# outfile = joinpath(outdir, "lapusta_kernel_fit_$(settings.nGroupsKII)_$(settings.nGroupsKIII)_$(settings.rhoWindow).jld2")
# JLD2.jldsave(outfile;
#     polesII=ComplexF64.(polesII),
#     residuesII=ComplexF64.(residuesII),
#     polesIII=ComplexF64.(polesIII),
#     residuesIII=ComplexF64.(residuesIII),
#     rhoWindow=Float64(settings.rhoWindow),
#     taperFraction=Float64(settings.taperFraction),
#     # Dense untapered exact table.  The brute-force stress law applies the
#     # same tail taper at setup for each q and fixed history time step.
#     exactRho=Float64.(exact.rho),
#     exactKII=Float64.(exact.KII),
#     exactKIII=Float64.(exact.KIII),
#     nu=Float64(settings.nu),
#     cp_cs=Float64(exact.cp_cs))
# @info "saved fit" outfile=outfile
