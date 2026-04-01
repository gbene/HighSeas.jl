using HighSeas
using CairoMakie
using GLMakie
using Peaks
using Glob
using CSV
using GLMakie.Colors
using Fractalizer
using Statistics
using StatsBase
using EasyFit
using CUDA

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


function GR(data)

    data_sort = sort(data,rev=true)
    N = log10.(1:length(data_sort))
    return data_sort, N/maximum(N)
end

input_dict = readSheet("../GPU_CUDA/BP4input.txt")
reduced_input_dict = readSheet("../GPU_CUDA/BP4input_reduced.txt", factor=14)
c_input_dict = copy(input_dict)
points = [[3e4, 3e4, -3e4, -3e4, 3e4] [1.5e4, -1.5e4, -1.5e4, 1.5e4, 1.5e4]] #RW patch points
np1 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[6442887322735277629, -8987213128142308954, -333252884332351366, 8464945807962482870])
np2 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[8513830690257299402, 5301462472252722888, 3588050925270478459, 3787793271851686014])
np3 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[6278251954719919174, 7730772275841400182, -2768929581059409545, 9035333999970238256])
np4 = NoiseParams(0.05:0.01:0.07, 1.0:1:10, -10.0:1:10.0, 100, 4, 10, seeds=[3947142509679233647, 626970696438225036, 7033574499115343085, -1445958421458504479])

templates = [random_template(n) for n in [np1,np2,np3,np4]]

shape = ClosedShape(points)
scaling = [1+2*input_dict["h"]/shape.l, 1+2*input_dict["h"]/shape.w]

buffer = shape * scaling

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
# samplers = loadData("temp/2026_03_02T13_10_19.489/CUDA/saved_SamplerSaver.jld2") # frac 4
# samplers = loadData("temp/2026_03_02T13_10_19.489/CUDA/saved_SamplerSaver.jld2") # frac 5
# samplers = loadData("../GPU_CUDA/BP4QD_out/2026_03_23T10_19_03.434/CUDA/saved_SamplerSaver.jld2")

# # samplers = loadData("../GPU_CUDA/BP4QD_out/2026_03_27T14_31_03.940/CUDA/saved_SamplerSaver.jld2")
# # samplers = loadData("../GPU_CUDA/BP4QD_out/2026_03_27T14_50_55.616/CUDA/saved_SamplerSaver.jld2")
# samplers = loadData("../GPU_CUDA/BP4QD_out/2026_03_27T15_00_26.360/CUDA/saved_SamplerSaver.jld2")

samplers = loadData("../GPU_CUDA/BP4QD_out/2026_03_31T16_06_16.866/CUDA/saved_SamplerSaver.jld2")


# catalog = loadData("../GPU_CUDA/BP4QD_out/2026_03_18T15_15_51.609/CUDA/saved_CatalogSaver.jld2")
# catalog = loadData("../GPU_CUDA/BP4QD_out/2026_03_27T14_31_03.940/CUDA/saved_CatalogSaver.jld2")
# catalog = loadData("../GPU_CUDA/BP4QD_out/2026_03_27T14_50_55.616/CUDA/saved_CatalogSaver.jld2")
catalog = loadData("../GPU_CUDA/BP4QD_out/2026_03_27T15_00_26.360/CUDA/saved_CatalogSaver.jld2")


bench8 = catalog
bench9 = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_02_19T12_50_10.812/CUDA/saved_CatalogSaver.jld2")
bench10 = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_02_19T14_41_08.325/CUDA/saved_CatalogSaver.jld2")
bench11 = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_02_21T15_29_51.441/CUDA/saved_CatalogSaver.jld2")
bench12 = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_02_23T15_40_29.769/CUDA/saved_CatalogSaver.jld2")
bench12f = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_fractal/2026_03_04T19_27_31.574/CUDA/saved_CatalogSaver.jld2")

bench12r_a = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_reduced_out/2026_03_20T12_47_55.139/CUDA/saved_CatalogSaver.jld2")
bench12r_b = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_reduced_out/2026_03_22T18_51_26.035/CUDA/saved_CatalogSaver.jld2")

bench12r = bench12r_a + bench12r_b

# bench13a = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026-02-11T15:07:57.681/CUDA/saved_CatalogSaver.jld2")
# bench13b = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026-02-12T11:28:13.033/CUDA/saved_CatalogSaver.jld2")
# bench13c = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026-02-14T15:22:07.937/CUDA/saved_CatalogSaver.jld2")
bench13d = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_02_26T15_40_22.736/CUDA/saved_CatalogSaver.jld2")
bench13e = loadSSH(url, username, private_file, public_file, "RateState/BP4QD/Julia/BP4QD_out/2026_03_05T15_46_48.027/CUDA/saved_CatalogSaver.jld2")




bench13 = bench13d+bench13e#bench13a+bench13b+bench13c+bench13d


event_dict = Dict(["0.04" => (500, bench8),"0.02" => (250, bench9), "0.01" => (125, bench10), "0.005" => (62.5, bench11), "0.0025" => (31.25, bench12), "0.0007"=>(12.5, bench13), "0.0001" => (1.5, bench12r)])

k_ordered = sort!(collect(keys(event_dict)), rev=true)
k = ["0.04", "0.005", "0.02", "0.0025", "0.01", "0.0007", "0.0001"]# Order is weird to have left to right reading order #sort!(collect(keys(event_dict)), rev=true)
real_idx = [1, 4, 2, 5, 3, 6, 7]


min_mag = zeros(6)

barbot_data = CSV.File(open("../resources/scec/barbot.6/b6_rupture.csv"))
cheng_data = CSV.File(open("../resources/scec/cheng/contour.csv"))

lambert_data = CSV.File(open("../resources/scec/lambert/contour.csv"))
ozawa_data = CSV.File(open("../resources/scec/ozawa/contour.csv"))


barbot_data.t[barbot_data.t .== 0.] .= 1e9

sampfig, sampax = HighSeas.plotPointSample(samplers.samplers[sample_point], "../resources/scec/", "slip_rate_2", "Vs", log10)

secfig, secax = HighSeas.plotSection(samplers.samplers[15], grid, 50, 150)

catalogfig, catalogax = HighSeas.plotCatalog(catalog, "Moment";ax_kwargs=(), stem_kwargs=())

domainfig, domainax = HighSeas.plotDomain(domain, samplers.samplers[sample_point], samplers.samplers[15])

areaslipfig = Figure(size=(1924,1080))
areaslipax  = Axis(areaslipfig[1,1], xlabel="Dc [m]", xticks=(1:length(k)-1, k_ordered[1:end-1]))
hideydecorations!(areaslipax, ticks = true)
hidexdecorations!(areaslipax, ticks = false, ticklabels=false, label=false, grid=true)



magfullfig = Figure(size=(1924,1080))
magfullax  = Axis(magfullfig[1,1], xlabel="Dc", ylabel="Magnitudes", xticks=(1:length(k)-1, k_ordered[1:end-1]))
magpartfig = Figure(size=(1924,1080))
magpartax  = Axis(magpartfig[1,1], xlabel="Dc", ylabel="Magnitudes", xticks=(1:length(k)-1, k_ordered[1:end-1]))
eventfig = Figure(size=(1924,1080))
slipfig = Figure(size=(1924,1080))

grfig = Figure(size=(1924, 1080))
grax = Axis(grfig[1,1], aspect=DataAspect(), title="GR distribution for partial ruptures", ylabel="Dc [m]", xlabel="Magnitudes", yticks=([1, 2.5, 4], ["0.005", "0.0025", "0.0007"]))

# hideydecorations!(grax)

indices = CartesianIndices(zeros(2,3))
areathresh = 0.8

bench12fig, bench12ax = HighSeas.plotCatalog(bench12,"mag";ax_kwargs=(), stem_kwargs=(markersize=0,stemcolor=Makie.wong_colors()[4]))
bench12rfig, bench12rax = HighSeas.plotCatalog(bench12r,"mag";ax_kwargs=(), stem_kwargs=(markersize=0,stemcolor=Makie.wong_colors()[7]))

for i in eachindex(k)
    dc = k[i]
    gs, b = event_dict[dc]

    if dc != "0.0001"

        c_input_dict["cellsizex"] = gs
        c_input_dict["cellsizey"] = gs
        localgrid = PowerGrid(c_input_dict)
        localfault = RectangleFault(c_input_dict, localgrid)

        localpatch = RectanglePatch(c_input_dict, localgrid)

        total_area = (sum(localpatch.dRW))*localgrid.cell_area
        # println("Area: $total_area")
        areas = b.Area/total_area
        times = b.t/(365*24*60*60)
        timemask = 158 .< times .< 652.0
        nanmask = .!isnan.(b.mag)

        events = b.mag[timemask]
        min_mag[real_idx[i]] = minimum(events)


        full_mask = areas .>= areathresh
        partial_mask = areas .< areathresh

        full_ruptures = areas[full_mask.*timemask] #full ruptures within the time mask
        full_slips = b.MeanSlip[full_mask.*timemask]


        fullevents = events[full_mask[timemask]]
        partialevents = events[partial_mask[timemask]]
        boxplot!(magfullax, fill(real_idx[i], length(fullevents)), fullevents, color=Makie.wong_colors()[i])
        boxplot!(magpartax, fill(real_idx[i], length(partialevents)), partialevents, color=Makie.wong_colors()[i])



        ax = Axis(eventfig[indices[i][1],indices[i][2]], title="Dc: $dc m")
        axslip = Axis(slipfig[indices[i][1],indices[i][2]], title="Dc: $dc m")

        if parse(Float64, dc) < 0.01
            mags, norm_vals = GR(partialevents)

            fig_mags = copy(mags)
            fit_norm_vals = copy(norm_vals)
            mask_mags = mags .>= minimum(mags)
            if parse(Float64, dc) == 0.0007
                mask_mags = 5.25 .< mags .< 6.35
            elseif parse(Float64, dc) == 0.0025
                mask_mags = 5.5 .< mags
            elseif parse(Float64, dc) == 0.005
                mask_mags = 6 .< mags
            end

            fit_norm_vals = norm_vals[mask_mags]
            fig_mags = mags[mask_mags]
            fit = fitlinear(fig_mags, fit_norm_vals)
            b_val = round(fit.a,digits=2)
            Rsq = round(fit.R2,digits=2)
            # display(fit)


            s = scatter!(grax, mags, norm_vals, color=Makie.wong_colors()[i])
            l = lines!(grax, fit.x, fit.y, color=Makie.wong_colors()[i], label="$dc", linestyle=:dash)
            a = annotation!(grax, fit.x[end], fit.y[end], fit.x[end], fit.y[end], text="b = $(-1*b_val), R² = $Rsq",align=(:left, :center))
            if parse(Float64, dc) == 0.0007
                mask_mags = 6.35 .< mags
                fit_norm_vals = norm_vals[mask_mags]
                fig_mags = mags[mask_mags]

                fit = fitlinear(fig_mags, fit_norm_vals)
                b_val = round(fit.a,digits=2)
                Rsq = round(fit.R2,digits=2)
                l2 = lines!(grax, fit.x, fit.y, color=Makie.wong_colors()[i], linestyle=:dash)
                a2 = annotation!(grax, fit.x[end], fit.y[end], fit.x[end], fit.y[end], text="b = $(-1*b_val), R² = $Rsq",align=(:left, :center))
                translate!(l2, 0, 3.5, 0)
                translate!(a2, 0, 3.5, 0)
            end


            # r = rainclouds!(grax, fill(dc, length(partialevents)), partialevents,markersize=5.0, orientation = :horizontal, color=Makie.wong_colors()[wang_color_idx[i]])
            # display(i)
            # translate!(s, 0, i/2, 0)
            # translate!(l, 0, i/2, 0)
            # translate!(a, 0, i/2, 0)
            if dc ==  "0.005"
                translate!(s, 0, 0.3, 0)
                translate!(l, 0, 0.5, 0)
                translate!(a, 0, 0.5, 0)
            elseif dc == "0.0025"
                translate!(s, 0, 1.3, 0)
                translate!(l, 0, 1.5, 0)
                translate!(a, 0, 1.5, 0)
            else
                translate!(s, 0, 3.3, 0)
                translate!(l, 0, 3.5, 0)
                translate!(a, 0, 3.5, 0)
            end


        end
        barplot!(areaslipax, [real_idx[i], real_idx[i]], [mean(full_ruptures), mean(full_slips)], dodge=[1, 2], bar_labels = :y, color=[Makie.wong_colors()[i], (Makie.wong_colors()[i],0.5)])


        HighSeas.plotCatalog(b, "mag", ax; stem_kwargs=(label="Dc: $dc",markersize=0, stemcolor=Makie.wong_colors()[i]))
        HighSeas.plotCatalog(b, "MeanSlip", axslip; stem_kwargs=(label="Dc: $dc",markersize=0, stemcolor=Makie.wong_colors()[i]))

        if i > 1

            HighSeas.plotCatalog(event_dict["0.04"][2], "mag", ax; stem_kwargs=(marker=:utriangle, color=Makie.wong_colors()[1],stemcolor=Makie.wong_colors()[1]))
            # HighSeas.plotCatalog(event_dict["0.04"][2], "MeanSlip", axslip; stem_kwargs=(marker=:utriangle, color=Makie.wong_colors()[1],stemcolor=Makie.wong_colors()[1]))

        end
        scatter!(ax, b.t[areas .>= areathresh]/(365*24*60*60), b.mag[areas .>= areathresh], marker=:utriangle, color=Makie.wong_colors()[i])
        scatter!(ax, b.t[areas .< areathresh]/(365*24*60*60), b.mag[areas .< areathresh], marker=:circle, color=Makie.wong_colors()[i])

        textlabel!(ax, (600, 4.2), text="N events: $(length(events))")

        scatter!(axslip, b.t[areas .>= areathresh]/(365*24*60*60), b.MeanSlip[areas .>= areathresh], marker=:utriangle, color=Makie.wong_colors()[i])
        scatter!(axslip, b.t[areas .< areathresh]/(365*24*60*60), b.MeanSlip[areas .< areathresh], marker=:circle, color=Makie.wong_colors()[i])



        if dc == "0.0025"

            HighSeas.plotCatalog(event_dict["0.04"][2], "mag", bench12ax; stem_kwargs=(marker=:utriangle, color=Makie.wong_colors()[1],stemcolor=Makie.wong_colors()[1]))

            scatter!(bench12ax, b.t[areas .>= areathresh]/(365*24*60*60), b.mag[areas .>= areathresh], marker=:utriangle, color=Makie.wong_colors()[i])
            scatter!(bench12ax, b.t[areas .< areathresh]/(365*24*60*60), b.mag[areas .< areathresh], marker=:circle, color=Makie.wong_colors()[i])

            textlabel!(bench12ax, (600, 3.2), text="N events: $(length(events))")

            fractal = fractalize(shape, templates, gs)

            localcpatch = CustomPatch(fractal.points, localgrid)

            cdomain = Domain(localgrid, localfault, localcpatch, nucleation)
            total_area = (sum(localcpatch.dRW))*localgrid.cell_area

            ylims!(bench12ax, 3, nothing)
            xlims!(bench12ax, 150, 660.0)

            # println("Fract area: $total_area")
            areas12f = bench12f.Area/total_area
            timesf = bench12f.t/(365*24*60*60)
            timemaskf = 158 .< timesf .< 652.0
            eventsf = bench12f.mag[timemaskf]
            fulleventsf = eventsf[areas12f[timemaskf] .>= areathresh]
            partialeventsf = eventsf[areas12f[timemaskf] .< areathresh]
            nanmaskf = .!isnan.(bench12f.mag)

            bench12fax = Axis(bench12fig[2,1])

            HighSeas.plotCatalog(bench12f,"mag",bench12fax;stem_kwargs=(markersize=0,stemcolor="#9fc4d8ff"))
            HighSeas.plotCatalog(event_dict["0.04"][2], "mag", bench12fax; stem_kwargs=(marker=:utriangle, color=Makie.wong_colors()[1],stemcolor=Makie.wong_colors()[1]))


            scatter!(bench12fax, bench12f.t[areas12f .>= areathresh]/(365*24*60*60), bench12f.mag[areas12f .>= areathresh], marker=:utriangle, color="#9fc4d8ff")
            scatter!(bench12fax, bench12f.t[areas12f .< areathresh]/(365*24*60*60), bench12f.mag[areas12f .< areathresh], marker=:circle, color="#9fc4d8ff")
            textlabel!(bench12fax, (600, 3.2), text="N events: $(length(eventsf))")

            ylims!(bench12fax, 3, nothing)
            xlims!(bench12fax, 150, 660.0)

            mags, norm_vals = GR(partialeventsf)

            fig_mags = copy(mags)
            fit_norm_vals = copy(norm_vals)
            mask_mags = 7 .> mags .>= 5

            fit_norm_vals = fit_norm_vals[mask_mags]
            fig_mags = mags[mask_mags]

            fit = fitlinear(fig_mags, fit_norm_vals)
            b_val = round(fit.a,digits=2)
            Rsq = round(fit.R2,digits=2)

            s = scatter!(grax, mags, norm_vals, color="#9fc4d8ff")
            l = lines!(grax, fit.x, fit.y, color="#9fc4d8ff", label="$dc Rough", linestyle=:dash)
            a = annotation!(grax, fit.x[end], fit.y[end], fit.x[end], fit.y[end], text="b = $(-1*b_val), R² = $Rsq",align=(:left, :center))

            translate!(s, 0, 2.3, 0)
            translate!(l, 0, 2.5, 0)
            translate!(a, 0, 2.5, 0)





            global bench12boxax = Axis(bench12fig[:,2], xlabel="Dc", ylabel="Magnitudes", xticks=([1,3], ["0.0025", "0.0007"]))


            boxplot!(bench12boxax, fill(0.5, length(partialevents)), partialevents,color=Makie.wong_colors()[i])
            boxplot!(bench12boxax, fill(1.5, length(partialeventsf)), partialeventsf, color="#9fc4d8ff")

        end

        if i == 6
            boxplot!(bench12boxax, fill(3, length(partialevents)), partialevents, color=Makie.wong_colors()[i])
        end
        println(i)
        ylims!(ax, 4, nothing)
        xlims!(ax, 150, 660.0)
        xlims!(axslip, 150, 660.0)
    else
        localgrid = PowerGrid(reduced_input_dict)
        localfault = RectangleFault(reduced_input_dict, localgrid)

        localpatch = RectanglePatch(reduced_input_dict, localgrid)

        total_area = (sum(localpatch.dRW))*localgrid.cell_area
        # println("Area: $total_area")
        areas = b.Area/total_area
        scatter!(bench12rax, b.t[areas .>= areathresh]/(365*24*60*60), b.mag[areas .>= areathresh], marker=:utriangle, color=Makie.wong_colors()[7])
        scatter!(bench12rax, b.t[areas .< areathresh]/(365*24*60*60), b.mag[areas .< areathresh], marker=:circle, color=Makie.wong_colors()[7])
        HighSeas.plotCatalog(event_dict["0.04"][2], "mag", bench12rax; stem_kwargs=(marker=:utriangle, color=Makie.wong_colors()[1],stemcolor=Makie.wong_colors()[1]))
        ylims!(bench12rax, 2, 5)
        # xlims!(bench12rax, 150, 660.0)

    end
end
println("asdasd")
axislegend(grax)
xlims!(grax, high=7.5)

fitfig = Figure(size=(1924,1080))
fitax = Axis(fitfig[1,1], xlabel="log10(Dc)", ylabel="Magnitudes")

dc_values = log10.(parse.(Float64, k_ordered[1:end-1]))
dc_pred = fitlinear(dc_values, min_mag)

Rsq_dc = round(dc_pred.R2,digits=2)
m_dc = round(dc_pred.a,digits=2)

extr = log10.(0.0001:0.0001:0.0007)

scatter!(fitax, dc_values, min_mag, label="Run models")

scatter!(fitax, log10(0.0001), minimum(bench12r.mag[bench12r.t/(360*24*60*60) .>20]), color=:red, marker=:utriangle, label="Reduced")
annotation!(fitax, dc_pred.x[end], dc_pred.y[end], dc_pred.x[end], dc_pred.y[end], text="m= $m_dc, R² = $Rsq_dc",align=(:left, :center))

lines!(fitax, dc_pred.x, dc_pred.y)
lines!(fitax, extr, dc_pred.(extr), linestyle=:dash)
axislegend(fitax)



contourfig, contourax = HighSeas.plotDomain(domain)
contour_data = Matrix(samplers.samplers[16].contour)
contour_data[isnan.(contour_data)] .= 1e9
contour!(contourax, grid.X, grid.Y, contour_data', labels=false,levels=30:10:60, color=:black, label="Ours")
contour!(contourax, barbot_data.x2, (barbot_data.x3.*-1), barbot_data.t, labels=false, levels=30:10:60, label="barbot.6", color=colorant"#0072b2ff")
contour!(contourax, cheng_data.x2, (cheng_data.x3.*-1), cheng_data.t, labels=false, levels=10:10:30, label="cheng", color=colorant"#e69f00ff")
contour!(contourax, lambert_data.x2, (lambert_data.x3.*-1).-500, lambert_data.t, labels=false, levels=30:10:60, label="lambert", color=colorant"#009e73ff" )
contour!(contourax, ozawa_data.x2, (ozawa_data.x3.*-1), ozawa_data.t, labels=false, levels=30:10:60, label="ozawa", color=colorant"#cc79a7ff")
xlims!(contourax, -35000, 35000)
ylims!(contourax, -20000, 20000)

interfig, interax = plotInterEventTime(samplers.samplers[sample_point], "../resources/scec/", "slip_rate_2", "Vs", log10)

peakfig, peakax = plotEventComparison(samplers.samplers[sample_point], "../resources/scec/", "slip_rate_2", "Vs", log10, ref_sim=3)
