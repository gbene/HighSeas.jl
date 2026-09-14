using GLMakie
using FFTW

ρw = 200
ρf = 0.2
ρs = (1-ρf)*ρw

ρ = 0:0.01:360



tapering_f(x, ρs, ρw) = 0.5 * (1 + cos(π * (x - ρs)/(ρw-ρs)))

tapering_f2(x, ρf, ρw) = 0.5 * (1 - cos(π * (x - ρw)/(ρw*ρf)))

tapering_hann(x, ρf, ρw) = 0.5 * (1 - cos(π * (ρw - x) / (ρw * ρf)))

# tapering_tukey(x, ρf, ρw, α) = 0 ≤ ρw - x ≤ ρw ? 0.5 * (1- cos(2π * (ρw - x) / (α * (ρw * ρf)))) : NaN

function calculate_response(window)

    # fft_window = zeros(2048)
    # fft_window[1:length(window)] = window

    A = fft(window)/(length(window)/2)

    freq = LinRange(-0.5, 0.5, length(A))

    response = log10.(abs.(fftshift( A / maximum(abs.(A)))))

    return freq, response
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

function rasied_cos(x, tapering_end, alpha)

    center = tapering_end/(1+alpha)

    tapering_start = center*(1-alpha)
    # tapering_end   =  center*(1+alpha)

    if x < tapering_start

        return 1

    elseif tapering_start ≤ x ≤ tapering_end

        if alpha == 0
            return 1
        else
            return 0.5 * ( 1 + cos( π * (x - tapering_start) / (2 * alpha * center) ) )
        end
    else

        return 0

    end
end




# fig, ax, plt = lines(ρ, tapering_f.(ρ, ρs, ρw))
# lines!(ax, ρ, tapering_f2.(ρ, ρf, ρw))
# lines!(ax, ρ, tapering_hann.(ρ, ρf, ρw))

alpha = 0.2
my_window = rasied_cos.(ρ, 200, alpha)
chat_window = tail_taper_window(ρ, 200, alpha)

my_freq, my_resp = calculate_response(my_window)
chat_freq, chat_resp = calculate_response(chat_window)

fig = Figure()

ax = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

lines!(ax2, my_freq, my_resp)
lines!(ax2, chat_freq, chat_resp)

lines!(ax, ρ, my_window)
lines!(ax, ρ, chat_window)

display(fig)
