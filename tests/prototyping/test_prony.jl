using Prony
using GLMakie
using SpecialFunctions
using Makie

x = 0:1e-1:10
y = @. besselj0(π*x)


fig, ax, plt = lines(x,y, label="data",linestyle=:solid)

pronyfunc = prony(x, y, 10)
xnew = LinRange(0,4,1000)
lines!(ax, xnew, real(pronyfunc.(xnew)), label="Prony fit",linestyle=:dash)

fig
