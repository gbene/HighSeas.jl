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




url = ENV["elja_url"]
username = ENV["elja_user"]
private_file = ENV["elja_private"]
public_file = ENV["elja_pub"]

# 256

# catalog = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_29T15_27_54.824/CUDA/saved_CatalogSaver.jld2")
# log_file = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_29T15_27_54.824/CUDA/simulation.log")

# 512

# catalog = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_29T15_33_31.501/CUDA/saved_CatalogSaver.jld2")
# log_file = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_29T15_33_31.501/CUDA/simulation.log")

# 1024

# catalog = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_29T16_23_59.628/CUDA/saved_CatalogSaver.jld2")
# log_file = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_29T16_23_59.628/CUDA/simulation.log")

# 2048

# catalog = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_29T17_22_23.321/CUDA/saved_CatalogSaver.jld2")
# log_file = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_29T17_22_23.321/CUDA/simulation.log")

# 4096

catalog1 = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_13T18_22_50.535/CUDA/saved_CatalogSaver.jld2")
catalog2 = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_20T21_40_15.145/CUDA/saved_CatalogSaver.jld2")
catalog3 = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_27T17_40_14.197/CUDA/saved_CatalogSaver.jld2")
catalog4 = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_29T23_13_54.261/CUDA/saved_CatalogSaver.jld2")
catalog5 = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_06_05T08_35_23.589/CUDA/saved_CatalogSaver.jld2")

catalog = catalog1 + catalog2 + catalog3 + catalog4 + catalog5
to_csv(catalog, "/home/gab28/DATA/PhD/GitHub/HighSeas.jl/examples/GPU_CUDA/Asperity_out/csvs/catalog_4096.csv")

log_file = loadSSH(url, username, private_file, public_file, "RateState/FractalAsperity/Asperity_out/2026_05_13T18_22_50.535/CUDA/simulation.log")

input_dict = log_file.input_dict
# Define the domain

grid = PowerGrid(input_dict)
fault = HighSeas.CircularFault(input_dict, grid, 0.8)
patch = HighSeas.FractalAsperity(input_dict, grid, 0.9, 0.35,16)

domain = Domain(grid, fault, patch)

figdomain, axdomain = plotDomain(domain, figdisplay=false)

figcatalog, axcatalog  = HighSeas.plotCatalog(catalog, "mag", figdisplay=true)

xlims!(axcatalog, 0, 3.1)
