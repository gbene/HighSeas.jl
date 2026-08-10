using HighSeas
using CairoMakie
using GLMakie
using Peaks
using Glob
using CSV
using GLMakie.Colors
using Statistics
using StatsBase
using CUDA


function GR(data)

    data_sort = sort(data,rev=true)
    N = log10.(1:length(data_sort))
    return data_sort, N/maximum(N)
end



# 128

catalog = loadData("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/Asperity_out/2026_06_08T09_20_43.863/CUDA/saved_CatalogSaver.jld2")
log_file = HighSeas.loadLog("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/Asperity_out/2026_06_08T09_20_43.863/CUDA/simulation.log")
to_csv(catalog, "/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/Asperity_out/2026_06_08T09_20_43.863/CUDA/catalog.csv")

input_dict = log_file.input_dict
# Define the domain

grid = PowerGrid(input_dict)
fault = HighSeas.CircularFault(input_dict, grid, 0.8)
patch = HighSeas.FractalAsperity(input_dict, grid, 0.9, 0.35,16)

domain = Domain(grid, fault, patch)

figdomain, axdomain = plotDomain(domain, figdisplay=false)

figcatalog, axcatalog  = HighSeas.plotCatalog(catalog, "mag", figdisplay=true)

xlims!(axcatalog, 0, 3.1)
