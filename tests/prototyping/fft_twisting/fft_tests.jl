using FFTW
using GLMakie


t = 1 #duration of signal in s

fs = 1024 # sampling frequency Hz = samp/s
T = 1/fs  # sampling period

N = fs*t # number of samples

t = range(0, t, N)

f = sin.(2π*t)
g = cos.(10*2π*t)

f_hat = fft(f) .* fft(g)

F = ifft(f_hat)

freq = fftfreq(N, fs)


fig = Figure()

ax = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
ax3 = Axis(fig[3,1])


lines!(ax, t, f)
lines!(ax2, fftshift(freq), abs.(fftshift(f_hat)))
lines!(ax3, t, real.(F))
display(fig)
