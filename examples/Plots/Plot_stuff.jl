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


input_dict = readSheet("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/circular_input.txt")
# saved_step = loadData("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/Asperity_out/2026_05_15T13_03_01.128/CUDA/saved_StepSaver.jld2")
# Define the domain

grid = PowerGrid(input_dict)
fault = HighSeas.CircularFault(input_dict, grid, 0.8)
patch = HighSeas.FractalAsperity(input_dict, grid, 0.9, 0.35,16)

domain = Domain(grid, fault, patch)

figdomain, axdomain = plotDomain(domain, figdisplay=false)

url = ENV["elja_url"]
username = ENV["elja_user"]
private_file = ENV["elja_private"]
public_file = ENV["elja_pub"]


# catalog = loadData("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/Asperity_out/2026_05_15T18_23_48.063/CUDA/saved_CatalogSaver.jld2") # 256
# catalog = loadData("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/Asperity_out/2026_05_18T08_50_31.472/CUDA/saved_CatalogSaver.jld2") # 512
# catalog = loadData("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/Asperity_out/2026_05_19T09_47_15.385/CUDA/saved_CatalogSaver.jld2") # 1024
catalog = loadData("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/Asperity_out/2026_05_21T08_02_53.495/CUDA/saved_CatalogSaver.jld2") # 2048


# catalog1 = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_13T18_22_50.535/CUDA/saved_CatalogSaver.jld2") # 4096
# catalog2 = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_20T21_40_15.145/CUDA/saved_CatalogSaver.jld2") # 4096
# catalog = catalog1 + catalog2

figcatalog, axcatalog  = HighSeas.plotCatalog(catalog, "mag", figdisplay=true)

xlims!(axcatalog, 0,3.5)
