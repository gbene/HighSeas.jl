using GLMakie
using SpecialFunctions

t = 0:0.1:100

b = -π/4

n = 2

a = 2

f(t,b) = cos(t+b)

f2(t,b) = (ℯ^(im * (b + t)) + ℯ^(-im * (b + t)))/2

g(t, b, n) = t^n*cos(t+b)

g2(t, b, n) = (factorial(n-1) * t^n) / 2 * (ℯ^(im * (b + t))+ ℯ^(-im * (b + t)))


k(t, b, a) = t^a * cos(t + b)
k2(t, b, a) = ((gamma(a+1,2.677) * t^a) / 2) * (ℯ^(im * (b + t))+ ℯ^(-im * (b + t)))

fig, ax, plt = lines(t, f.(t, b))

lines!(ax, t, real(f2.(t, b)), linestyle=:dash)

fig

fig2, ax2, plt2 = lines(t, g.(t, b, n))

lines!(ax2, t, real(g2.(t, b, n)), linestyle=:dash)

fig2

fig3, ax3, plt3 = lines(t, k.(t, b, a))

lines!(ax3, t, real(k2.(t, b, a)), linestyle=:dash)

fig3
