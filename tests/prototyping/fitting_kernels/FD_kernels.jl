using GLMakie
using FromFile

@from "Lapustakernels.jl" import LapustaKernels.lapustakernels





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
