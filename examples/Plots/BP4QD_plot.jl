using HighSeas
using CairoMakie
using GLMakie
using JLD2
using Peaks
using Glob
using CSV
using GLMakie.Colors
using Fractalizer
using Statistics

function plotInterEventTime(pointSampler, ref_path::String, quantity::String, sampler_quantity::String, scale::Function; display=false)

    function calcInter(property, time)
        peaks = findmaxima(property,30).indices
        peaktime = time[peaks]/(365*24*60*60)
        interevent = peaktime[2:end]-peaktime[1:end-1]
        # pushfirst!(interevent, 0.0)
        return interevent
    end

    point = pointSampler.sample_point_id
    pointx = pointSampler.sample_point_x
    pointy = pointSampler.sample_point_y


    paths = glob("**/$point.csv", ref_path)
    data = Array{CSV.File, 1}(undef, length(paths))
    label = Array{String, 1}(undef, length(paths))
    fig = Figure(size=(1920,1080), figure_padding=30)


    ax = Axis(fig[1,1], title="Interevent time for point $point (x:$pointx, y:$pointy)", xlabel="Event n", ylabel="Interevent time [yr]")


    for i in eachindex(paths)
        path = paths[i]
        data = CSV.File(open(path))
        label = splitpath(path)[end-1]
        property = getproperty(data, Symbol(quantity))
        interevent = calcInter(property, data.t)
        lines!(ax, interevent, label=label, linewidth=5)
        scatter!(ax, interevent, label=label, markersize=15)

    end

    property = scale.(getproperty(pointSampler, Symbol(sampler_quantity)))
    interevent = calcInter(property, pointSampler.times)
    lines!(ax, interevent, label="Ours", color=:black, linewidth=5)
    scatter!(ax, interevent, label="Ours",color=:black, markersize=15)

    ax.titlesize=40
    ax.xlabelsize = 40
    ax.ylabelsize = 40
    ax.xticklabelsize = 40
    ax.yticklabelsize = 40
    axislegend(labelsize=40, merge=true, position=:rb)


    if display
        display(fig)
    end

    return fig, ax

end


function plotEventComparison(pointSampler, ref_path::String, quantity::String, sampler_quantity::String, scale::Function; display=false, ref_sim::Int=1)

    function calcPeakTime(property, time)
        peaks = findmaxima(property,30).indices
        peaktime = time[peaks]/(365*24*60*60)
        return peaktime
    end

    point  = pointSampler.sample_point_id
    pointx = pointSampler.sample_point_x
    pointy = pointSampler.sample_point_y


    paths = glob("**/$point.csv", ref_path)
    data = Array{CSV.File, 1}(undef, length(paths))
    label = Array{String, 1}(undef, length(paths))
    fig = Figure(size=(1920,1080), figure_padding=30)


    ax = Axis(fig[1,1], title="Peak time difference for point $point (x:$pointx, y:$pointy)", xlabel="Event n", ylabel="Peak time difference [yr]")

    peaktimes = zeros(length(paths)+1, 10)
    diff_mat = zeros(length(paths)+1, 10)
    for i in eachindex(paths)
        path = paths[i]
        data = CSV.File(open(path))
        label = splitpath(path)[end-1]
        property = getproperty(data, Symbol(quantity))
        peaktime = calcPeakTime(property, data.t)
        peaktimes[i, 1:length(peaktime)] = peaktime
    end

    property = scale.(getproperty(pointSampler, Symbol(sampler_quantity)))
    peaktime = calcPeakTime(property, pointSampler.times)
    peaktimes[end, 1:length(peaktime)] = peaktime

    for i in 1:size(peaktimes,1)
        ref = peaktimes[ref_sim, :]
        row = peaktimes[i, :]
        # println(i)
        if i == 5
            lines!(ax, row-ref, linewidth=5, label="Ours", color=:black)
            scatter!(ax, row-ref, markersize=15, label="Ours", color=:black)


        else
            path = paths[i]
            label = splitpath(path)[end-1]
            scatter!(ax, row-ref, markersize=15, label=label)
            lines!(ax, row-ref, linewidth=5, label=label)
        end
    end


    # lines!(ax, interevent, label="Ours", color=:black, linewidth=5)
    # scatter!(ax, interevent, label="Ours",color=:black, markersize=15)

    ax.titlesize=40
    ax.xlabelsize = 40
    ax.ylabelsize = 40
    ax.xticklabelsize = 40
    ax.yticklabelsize = 40
    ylims!(ax, -20,20)
    axislegend(labelsize=40, merge=true, position=:lb)


    if display
        display(fig)
    end

    return fig, ax

end


function plotContour(contourSampler, ref_path::String)


    paths = glob("**/contour.csv", ref_path)
    data = Array{CSV.File, 1}(undef, length(paths))
    label = Array{String, 1}(undef, length(paths))
    fig = Figure(size=(1920,1080), figure_padding=30)

    ax = Axis(fig[1,1], title="Firt rupture contour time", xlabel="X [m]", ylabel="Y[m]")


    for i in eachindex(paths)
        path = paths[i]
        data = CSV.File(open(path))
        label = splitpath(path)[end-1]
        property = getproperty(data, Symbol(quantity))
        interevent = calcInter(property, data.t)
        lines!(ax, interevent, label=label, linewidth=5)
        scatter!(ax, interevent, label=label, markersize=15)

    end

    property = scale.(getproperty(pointSampler, Symbol(sampler_quantity)))
    interevent = calcInter(property, pointSampler.times)
    lines!(ax, interevent, label="Ours", color=:black, linewidth=5)
    scatter!(ax, interevent, label="Ours",color=:black, markersize=15)

    ax.titlesize=40
    ax.xlabelsize = 40
    ax.ylabelsize = 40
    ax.xticklabelsize = 40
    ax.yticklabelsize = 40
    axislegend(labelsize=40, merge=true, position=:rb)


    if display
        display(fig)
    end

    return fig, ax
end


function GR(magnitudes, data)
    N = zeros(length(magnitudes))

    for i in eachindex(N)
        N[i] = sum(data.>=magnitudes[i])
    end
    return N
end

input_dict = ReadSheet("BP4input.txt")
c_input_dict = copy(input_dict)
points = [[3e4, 3e4, -3e4, -3e4, 3e4] [1.5e4, -1.5e4, -1.5e4, 1.5e4, 1.5e4]] #RW patch points
np1 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[6442887322735277629, -8987213128142308954, -333252884332351366, 8464945807962482870])
np2 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[8513830690257299402, 5301462472252722888, 3588050925270478459, 3787793271851686014])
np3 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[6278251954719919174, 7730772275841400182, -2768929581059409545, 9035333999970238256])
np4 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[3947142509679233647, 626970696438225036, 7033574499115343085, -1445958421458504479])

templates = [random_template(n) for n in [np1,np2,np3,np4]]

shape = ClosedShape(points)
s = [1+2*input_dict["h"]/shape.l, 1+2*input_dict["h"]/shape.w]

buffer = shape * s

url = ENV["elja_url"]
username = ENV["elja_user"]
private_file = ENV["elja_private"]
public_file = ENV["elja_pub"]

grid = PowerGrid(input_dict)
fault = RectangleFault(input_dict, grid)



patch = RectanglePatch(input_dict, grid)
nucleation = RectangleNucleation(input_dict, grid)

domain = Domain(grid, fault, patch, nucleation)

sample_point = 8
# samplers = loadData("temp/2026_03_02T13_19_10.357/CUDA/saved_SamplerSaver.jld2") # frac 2
# samplers = loadData("temp/2026_03_02T13_15_40.517/CUDA/saved_SamplerSaver.jld2") # frac 3
samplers = loadData("temp/2026_03_02T13_10_19.489/CUDA/saved_SamplerSaver.jld2") # frac 4
# samplers = loadData("temp/2026_03_02T13_10_19.489/CUDA/saved_SamplerSaver.jld2") # frac 5


catalog = loadData("BP4QD_out/2026_02_13T11_40_14.783/CUDA/saved_CatalogSaver.jld2")


bench8 = loadData("BP4QD_out/2026_02_24T11_36_34.660/CUDA/saved_CatalogSaver.jld2")
bench9 = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_02_19T12_50_10.812/CUDA/saved_CatalogSaver.jld2")
bench10 = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_02_19T14_41_08.325/CUDA/saved_CatalogSaver.jld2")
bench11 = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_02_21T15_29_51.441/CUDA/saved_CatalogSaver.jld2")
bench12 = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_02_23T15_40_29.769/CUDA/saved_CatalogSaver.jld2")
bench12f = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_fractal/2026_03_04T19_27_31.574/CUDA/saved_CatalogSaver.jld2")

# bench13a = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026-02-11T15:07:57.681/CUDA/saved_CatalogSaver.jld2")
# bench13b = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026-02-12T11:28:13.033/CUDA/saved_CatalogSaver.jld2")
# bench13c = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026-02-14T15:22:07.937/CUDA/saved_CatalogSaver.jld2")
bench13d = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_02_26T15_40_22.736/CUDA/saved_CatalogSaver.jld2")
bench13e = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_03_05T15_46_48.027/CUDA/saved_CatalogSaver.jld2")




bench13 = bench13d+bench13e#bench13a+bench13b+bench13c+bench13d


event_dict = Dict(["0.04" => (500, bench8),"0.02" => (250, bench9), "0.01" => (125, bench10), "0.005" => (62.5, bench11), "0.0025" => (31.25, bench12), "0.0007"=>(12.5, bench13)])

k = sort!(collect(keys(event_dict)), rev=true)


barbot_data = CSV.File(open("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/resources/scec/barbot.6/b6_rupture.csv"))
cheng_data = CSV.File(open("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/resources/scec/cheng/contour.csv"))

lambert_data = CSV.File(open("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/resources/scec/lambert/contour.csv"))
ozawa_data = CSV.File(open("/home/gab28/DATA/PhD/GitHub/HighSeas.jl/resources/scec/ozawa/contour.csv"))


barbot_data.t[barbot_data.t .== 0.] .= 1e9

sampfig, sampax = HighSeas.plotPointSample(samplers.samplers[sample_point], "../resources/scec/", "slip_rate_2", "Vs", log10)

secfig, secax = HighSeas.plotSection(samplers.samplers[15], 50, 150)

catalogfig, catalogax = HighSeas.plotCatalog(catalog, "Moment";ax_kwargs=(), stem_kwargs=())

domainfig, domainax = HighSeas.plotDomain(domain, samplers.samplers[sample_point], samplers.samplers[15])

areaslipfig = Figure(size=(1924,1080))
areaslipax  = Axis(areaslipfig[1,1], xlabel="Dc", xticks=(1:length(k), k[1:end]))
hideydecorations!(areaslipax, ticks = true)
hidexdecorations!(areaslipax, ticks = false, ticklabels=false, label=false, grid=true)



magfullfig = Figure(size=(1924,1080))
magfullax  = Axis(magfullfig[1,1], xlabel="Dc", ylabel="Magnitudes", xticks=(1:length(k), k[1:end]))
magpartfig = Figure(size=(1924,1080))
magpartax  = Axis(magpartfig[1,1], xlabel="Dc", ylabel="Magnitudes", xticks=(1:length(k), k[1:end]))
eventfig = Figure(size=(1924,1080))


indices = CartesianIndices(zeros(2,3))
areathresh = 0.9

bench12fig, bench12ax = HighSeas.plotCatalog(bench12,"mag";ax_kwargs=(), stem_kwargs=(markersize=0,))



for i in eachindex(k)

    dc = k[i]
    gs, b = event_dict[dc]

    times = b.t/(365*24*60*60)

    timemask = 158 .< times .< 652.0
    nanmask = .!isnan.(b.mag)

    events = b.mag[timemask]


    c_input_dict["cellsizex"] = gs
    c_input_dict["cellsizey"] = gs
    localgrid = PowerGrid(c_input_dict)
    localfault = RectangleFault(c_input_dict, localgrid)

    localpatch = RectanglePatch(c_input_dict, localgrid)

    total_area = (sum(localpatch.dRW)+sum(localpatch.dTR))*localgrid.cell_area

    areas = b.Area/total_area
    full_mask = areas .>= areathresh
    partial_mask = areas .< areathresh

    full_ruptures = areas[full_mask.*timemask]
    full_slips = b.MeanSlip[full_mask.*timemask]


    fullevents = events[full_mask[timemask]]
    partialevents = events[partial_mask[timemask]]
    boxplot!(magfullax, fill(i, length(fullevents)), fullevents)
    boxplot!(magpartax, fill(i, length(partialevents)), partialevents)



    ax = Axis(eventfig[indices[i][1],indices[i][2]], title="Dc: $dc")

    barplot!(areaslipax, [i, i], [mean(full_ruptures), mean(full_slips)], dodge=[1, 2], bar_labels = :y)


    HighSeas.plotCatalog(b, "mag", ax; stem_kwargs=(label="Dc: $dc",markersize=0))
    if i > 1

        HighSeas.plotCatalog(event_dict["0.04"][2], "mag", ax; stem_kwargs=(marker=:utriangle, color=(:orange, 0.5),stemcolor=(:orange, 0.5)))
    end
    scatter!(ax, b.t[areas .>= areathresh]/(365*24*60*60), b.mag[areas .>= areathresh], marker=:utriangle, color=Makie.wong_colors()[1])
    scatter!(ax, b.t[areas .< areathresh]/(365*24*60*60), b.mag[areas .< areathresh], marker=:circle, color=Makie.wong_colors()[1])

    if i == 5

        HighSeas.plotCatalog(event_dict["0.04"][2], "mag", bench12ax; stem_kwargs=(marker=:utriangle, color=(:orange, 0.5),stemcolor=(:orange, 0.5)))

        scatter!(bench12ax, b.t[areas .>= areathresh]/(365*24*60*60), b.mag[areas .>= areathresh], marker=:utriangle, color=Makie.wong_colors()[1])
        scatter!(bench12ax, b.t[areas .< areathresh]/(365*24*60*60), b.mag[areas .< areathresh], marker=:circle, color=Makie.wong_colors()[1])

        bench12fax = Axis(bench12fig[2,1])

        HighSeas.plotCatalog(bench12f,"mag",bench12fax;stem_kwargs=(markersize=0,))
        HighSeas.plotCatalog(event_dict["0.04"][2], "mag", bench12fax; stem_kwargs=(marker=:utriangle, color=(:orange, 0.5),stemcolor=(:orange, 0.5)))

        fractal = fractalize(shape, templates, gs)

        localcpatch = CustomPatch(fractal.points, localgrid)

        cdomain = Domain(localgrid, localfault, localcpatch, nucleation)
        total_area = (sum(localcpatch.dRW))*localgrid.cell_area
        areas12f = bench12f.Area/total_area

        scatter!(bench12fax, bench12f.t[areas12f .>= areathresh]/(365*24*60*60), bench12f.mag[areas12f .>= areathresh], marker=:utriangle, color=Makie.wong_colors()[1])
        scatter!(bench12fax, bench12f.t[areas12f .< areathresh]/(365*24*60*60), bench12f.mag[areas12f .< areathresh], marker=:circle, color=Makie.wong_colors()[1])

        ylims!(bench12fax, 3, nothing)
        xlims!(bench12fax, 150, 660.0)
        ylims!(bench12ax, 3, nothing)
        xlims!(bench12ax, 150, 660.0)

        eventsf = bench12f.mag[timemask]
        fulleventsf = eventsf[areas12f[timemask] .>= areathresh]
        partialeventsf = eventsf[areas12f[timemask] .< areathresh]


        global bench12boxax = Axis(bench12fig[:,2], xlabel="Dc", ylabel="Magnitudes", xticks=(1:3, ["normal", "fractal", "0.0007"]))


        boxplot!(bench12boxax, fill(1, length(partialevents)), partialevents)
        boxplot!(bench12boxax, fill(2, length(partialeventsf)), partialeventsf)

    end

    if i == 6
        boxplot!(bench12boxax, fill(3, length(partialevents)), partialevents)
    end

    ylims!(ax, 4, nothing)
    xlims!(ax, 150, 660.0)
end












contourfig, contourax = HighSeas.plotDomain(domain)
contour_data = Matrix(samplers.samplers[16].contour)
contour_data[isnan.(contour_data)] .= 1e9
contour!(contourax, grid.X, grid.Y, contour_data', labels=true,levels=30:10:60, color=:black, label="Ours")
contour!(contourax, barbot_data.x2, (barbot_data.x3.*-1), barbot_data.t, labels=false, levels=30:10:60, label="barbot.6", color=colorant"#0072b2ff")
contour!(contourax, cheng_data.x2, (cheng_data.x3.*-1), cheng_data.t, labels=false, levels=10:10:30, label="cheng", color=colorant"#e69f00ff")
contour!(contourax, lambert_data.x2, (lambert_data.x3.*-1).-500, lambert_data.t, labels=false, levels=30:10:60, label="lambert", color=colorant"#009e73ff" )
contour!(contourax, ozawa_data.x2, (ozawa_data.x3.*-1), ozawa_data.t, labels=false, levels=30:10:60, label="ozawa", color=colorant"#cc79a7ff")
xlims!(contourax, -35000, 35000)
ylims!(contourax, -20000, 20000)

interfig, interax = plotInterEventTime(samplers.samplers[sample_point], "../resources/scec/", "slip_rate_2", "Vs", log10)

peakfig, peakax = plotEventComparison(samplers.samplers[sample_point], "../resources/scec/", "slip_rate_2", "Vs", log10, ref_sim=3)
