"""

Understanding the wieght function of the kernel fit

"""


using GLMakie
using FromFile
@from "Lapustakernels.jl" import LapustaKernels.lapustakernels


rho, K2, K2int, K3, K3int = lapustakernels(0.25, 180.0, 1e-3)
