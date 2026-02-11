using HighSeas
using CairoMakie
using GLMakie
using JLD2


input_dict = ReadSheet("BP4input.txt")

grid = PowerGrid(input_dict)
fault = RectangleFault(input_dict, grid)

patch = RectanglePatch(input_dict, grid)
nucleation = RectangleNucleation(input_dict, grid)

domain = Domain(grid, fault, patch, nucleation)

sample_point = 8
samplers = load("BP4QD_out/2026-02-11T17:05:54.264/CPU/saved_SamplerSaver.jld2")["data"]
catalog = load("BP4QD_out/2026-02-11T17:05:54.264/CPU/saved_CatalogSaver.jld2")["data"]
sampfig, sampax = HighSeas.plotPointSample(samplers[sample_point], "../resources/scec/", "slip_2", "dxs")

secfig, secax = HighSeas.plotSection(samplers[15], 50, 150)

catalogfig, catalogax = HighSeas.plotCatalog(catalog, "Moment";ax_kwargs=(), stem_kwargs=())

domainfig, domainax = HighSeas.plotDomain(domain, samplers[sample_point], samplers[15])
