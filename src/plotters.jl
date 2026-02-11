abstract type AbstractPlotter end
abstract type LivePlotter <: AbstractPlotter end

struct RSPlotter <: LivePlotter

    state::AbstractState
    stepper::AbstractStepper
    plot_every::Int
    fig::Figure
    axs::Vector{Axis}
    obs::Vector{Observable}
    label::Label

    function RSPlotter(experiment::AbstractExperiment, algorithm::AbstractAlgorithm, plot_every::Int)
        fig = Figure()
        ax1  = Axis(fig[1,1], title="Slip", aspect=DataAspect())
        ax2  = Axis(fig[1,3], title="Log10(τ)", aspect=DataAspect())
        ax3  = Axis(fig[2,1], title="Log10(V)", aspect=DataAspect())
        ax4  = Axis(fig[2,3], title="Log10(θ)", aspect=DataAspect())


        state = experiment.state
        stepper = algorithm.stepper
        x = experiment.domain.grid.X
        y = experiment.domain.grid.Y

        dx = state.dx
        tau = state.tau
        V = state.V
        theta = state.theta

        # This is very dumb. To display correctly the data we need to do the adjoint
        # If we do Observable(dx') then it defines an observable that looks at the adjoint
        # However if we then do log10.() then the output will be of type Matrix so we observe only dx as an adjoint
        # This however does not trigger the type warning for some reason..
        # Moreover when we update the observable it does not care if it is a matrix or adjoint it just considers them as martrices thus if
        # we don't do the following the plotting will not work

        obs1 = Observable(dx)
        obs1[] = dx'
        obs2 = Observable(tau)
        obs2[] = log10.(tau')
        obs3 = Observable(V)
        obs3[] = log10.(V')
        obs4 = Observable(theta)
        obs4[] = log10.(theta')
        L = Label(fig[0,:], "Simulation time: 0.0sec\n Step: 0", fontsize=15)

        axs = [ax1,ax2,ax3,ax4]
        obs = [obs1,obs2,obs3,obs4]

        hm1 = heatmap!(ax1, x, y, obs1, rasterize=10)
        hm2 = heatmap!(ax2, x, y, obs2, rasterize=10)
        hm3 = heatmap!(ax3, x, y, obs3, rasterize=10)
        hm4 = heatmap!(ax4, x, y, obs4, rasterize=10)


        Colorbar(fig[1,2], hm1)
        Colorbar(fig[1,4], hm2)
        Colorbar(fig[2,2], hm3)
        Colorbar(fig[2,4], hm4)

        display(fig)


        new(state, stepper, plot_every, fig, axs, obs, L)
    end


end


function UpdatePlot(rsPlotter)
    step = rsPlotter.stepper.step

    dx = rsPlotter.state.dx
    tau = rsPlotter.state.tau
    V = rsPlotter.state.V
    theta = rsPlotter.state.theta

    if step % rsPlotter.plot_every == 0

        rsPlotter.obs[1][] = dx'
        rsPlotter.obs[2][] = log10.(tau')
        rsPlotter.obs[3][] = log10.(V')
        rsPlotter.obs[4][] = log10.(theta)'
        rsPlotter.label.text = "Simulation time: $(round(rsPlotter.stepper.time, digits=2))sec\n Step: $step"

        yield()
    end

end


function plotPointSample(pointSampler::PointSampler, sampler_quantity::String, display=false)

    point = pointSampler.sample_point_id
    pointx = pointSampler.sample_point_x
    pointy = pointSampler.sample_point_y

    fig = Figure(size=(1920,1080))


    ax = Axis(fig[1,1], title="Sample point $point comparison (x:$pointx, y:$pointy)", xlabel="Time [yr]", ylabel="$sampler_quantity")


    lines!(ax, pointSampler.times/(365*24*60*60), getproperty(pointSampler, Symbol(sampler_quantity)), label="Ours", color=:black)

    ax.titlesize=30
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.xticklabelsize = 25
    ax.yticklabelsize = 25
    axislegend()

    if display
        display(fig)
    end

    return fig, ax

end
function plotPointSample(pointSampler::PointSampler, sampler_quantity::String, scale::Function, display=false)

    point = pointSampler.sample_point_id
    pointx = pointSampler.sample_point_x
    pointy = pointSampler.sample_point_y

    fig = Figure(size=(1920,1080))


    ax = Axis(fig[1,1], title="Sample point $point comparison (x:$pointx, y:$pointy)", xlabel="Time [yr]", ylabel="$(string(scale))($sampler_quantity)")


    lines!(ax, pointSampler.times/(365*24*60*60), scale.(getproperty(pointSampler, Symbol(sampler_quantity))), label="Ours", color=:black)

    ax.titlesize=30
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.xticklabelsize = 25
    ax.yticklabelsize = 25
    axislegend()

    if display
        display(fig)
    end

    return fig, ax

end
function plotPointSample(pointSampler::PointSampler, ref_path::String, quantity::String, sampler_quantity::String, display=false)

    point = pointSampler.sample_point_id
    pointx = pointSampler.sample_point_x
    pointy = pointSampler.sample_point_y


    paths = glob("**/$point.csv", ref_path)
    data = Array{CSV.File, 1}(undef, length(paths))
    label = Array{String, 1}(undef, length(paths))
    fig = Figure(size=(1920,1080), figure_padding=30)


    ax = Axis(fig[1,1], title="Sample point $point comparison (x:$pointx, y:$pointy)", xlabel="Time [yr]", ylabel="$sampler_quantity")


    for i in eachindex(paths)
        path = paths[i]
        data = CSV.File(open(path))
        label = splitpath(path)[end-1]
        lines!(ax, data.t/(365*24*60*60), getproperty(data, Symbol(quantity)), label=label, linewidth=5)
    end

    lines!(ax, pointSampler.times/(365*24*60*60), getproperty(pointSampler, Symbol(sampler_quantity)), label="Ours", color=:black, linewidth=5)

    ax.titlesize=40
    ax.xlabelsize = 40
    ax.ylabelsize = 40
    ax.xticklabelsize = 40
    ax.yticklabelsize = 40
    axislegend(labelsize=40)


    if display
        display(fig)
    end
    return fig, ax

end
function plotPointSample(pointSampler::PointSampler, ref_path::String, quantity::String, sampler_quantity::String, scale::Function, display=false)

    point = pointSampler.sample_point_id
    pointx = pointSampler.sample_point_x
    pointy = pointSampler.sample_point_y


    paths = glob("**/$point.csv", ref_path)
    data = Array{CSV.File, 1}(undef, length(paths))
    label = Array{String, 1}(undef, length(paths))
    fig = Figure(size=(1920,1080), figure_padding=30)


    ax = Axis(fig[1,1], title="Sample point $point comparison (x:$pointx, y:$pointy)", xlabel="Time [yr]", ylabel="$(string(scale))($sampler_quantity)")


    for i in eachindex(paths)
        path = paths[i]
        data = CSV.File(open(path))
        label = splitpath(path)[end-1]
        lines!(ax, data.t/(365*24*60*60), getproperty(data, Symbol(quantity)), label=label, linewidth=5)
    end

    lines!(ax, pointSampler.times/(365*24*60*60), scale.(getproperty(pointSampler, Symbol(sampler_quantity))), label="Ours", color=:black, linewidth=5)

    ax.titlesize=40
    ax.xlabelsize = 40
    ax.ylabelsize = 40
    ax.xticklabelsize = 40
    ax.yticklabelsize = 40
    axislegend(labelsize=40)


    if display
        display(fig)
    end

    return fig, ax

end

function plotSection(sectionSampler::SectionSampler, s_slip_freq=50, as_slip_freq=50, display=false)


    X = sectionSampler.grid.X
    Y = sectionSampler.grid.Y
    section = sectionSampler.section
    seismic_slips = section[section[:,1] .> 1e-2, 2:end]
    aseismic_slips = section[section[:,1] .<= 1e-7, 2:end]

    fig = Figure(size=(480,1080))
    if sectionSampler.axis == "x"
        ax = Axis(fig[1,1], title="Along dip slip profile (x = $(sectionSampler.coord))",xlabel="Dip distance [Km]", ylabel="Cumulative slip")
        series!(ax, Y/1000, seismic_slips[1:s_slip_freq:end,:], solid_color=:red, linewidth=1)
        series!(ax, Y/1000, aseismic_slips[1:as_slip_freq:end,:], solid_color=:blue, linewidth=1)
    else
        ax = Axis(fig[1,1], title="Along strike slip profile (y = $(sectionSampler.coord))",xlabel="Strike distance [Km]", ylabel="Cumulative slip")
        series!(ax, X/1000, seismic_slips[1:s_slip_freq:end,:], solid_color=:red, linewidth=1)
        series!(ax, X/1000, aseismic_slips[1:as_slip_freq:end,:], solid_color=:blue, linewidth=1)
    end

    tightlimits!(ax)
    resize_to_layout!(fig)

    ax.titlesize=20
    ax.xlabelsize = 20
    ax.ylabelsize = 20
    ax.xticklabelsize = 20
    ax.yticklabelsize = 20

    legend_seismic = LineElement(color=:red, linewidth=1)
    legend_aseismic = LineElement(color=:blue, linewidth=1)

    axislegend(ax, [legend_seismic, legend_aseismic], ["Seismic", "Aseismic"])

    if display
        display(fig)
    end

    return fig, ax
end


function plotDomain(domain::AbstractDomain, dot_grid=false, display=false)
    fig = Figure(size=(480,480))
    ax = Axis(fig[1,1], aspect = DataAspect(), xlabel= "X [m]", ylabel= "Y [m]", title="Simulated domain")

    grid  = domain.grid
    fault = domain.fault
    patch = domain.patch
    nucleation = domain.nucleation


    heatmap!(ax, grid.X, grid.Y, Matrix(patch.dRS'),colormap=[(:black), :cyan], rasterize=10)

    heatmap!(ax, grid.X, grid.Y, Matrix(fault.dCR'),colormap=[(:white, 0.), :gray96], rasterize=10)
    # heatmap!(ax, grid.X, grid.Y, fault.dLO',colormap=[(:white, 0.), :cyan])


    heatmap!(ax, grid.X, grid.Y, Matrix(patch.dRW'),colormap=[(:white,0.), :springgreen1], rasterize=10)

    if patch.h > 0.0
        heatmap!(ax, grid.X, grid.Y, Matrix(patch.dTR'),colormap=[(:white,0.), :yellow], rasterize=10)
    end

    if ~(typeof(nucleation) <: EmptyNucleation)
        heatmap!(ax, grid.X, grid.Y, Matrix(nucleation.dNU'),colormap=[(:white,0.), :darkgreen], rasterize=10)
        # scatter!(ax, nucleation.xi, nucleation.yi, color=:black)
    end

    if dot_grid
        scatter!(ax, [(x, y) for x in grid.X for y in grid.Y], markersize=3,strokecolor=:white,strokewidth=0.5)
    end

    tightlimits!(ax)
    resize_to_layout!(fig)

    if display
        display(fig)
    end

    return fig, ax

end

function plotDomain(domain::AbstractDomain, sample_point::PointSampler, sample_section::SectionSampler, dot_grid=false, display=false; sample_point_kwargs=(), sample_section_kwargs=())

    fig, ax = plotDomain(domain, dot_grid)


    scatter!(ax, sample_point.sample_point_x, sample_point.sample_point_y; sample_point_kwargs...)

    if sample_section.axis == "y"
        hlines!(ax, sample_section.coord; sample_section_kwargs...)
    else
        vlines!(ax, sample_section.coord; sample_section_kwargs...)
    end

    if display
        display(fig)
    end
    return fig, ax

end

function plotCatalog(catalog::AbstractCatalog, quantity::String, display=false; ax_kwargs, stem_kwargs)

    n_events = catalog.n_events
    trimmed_catalog = Catalog(catalog.catalog[1:n_events, :])
    data = getfield(trimmed_catalog, Symbol(quantity))

    fig = Figure(size=(1920,1080), figure_padding=30)


    ax = Axis(fig[1,1], title="Simulated catalog", xlabel="Time [yrs]", ylabel=quantity; ax_kwargs...)

    stem!(ax, trimmed_catalog.t/(365*24*60*60), data; stem_kwargs...)

    ax.titlesize=30
    ax.xticklabelsize = 25
    ax.xlabelsize = 25
    ax.ylabelsize = 25
    ax.yticklabelsize = 25

    # axislegend(ax)

    if display
        display(fig)
    end

    return fig, ax
end
function plotCatalog(catalog::AbstractCatalog, quantity::String, ax::Axis; stem_kwargs)

    n_events = catalog.n_events
    trimmed_catalog = Catalog(catalog.catalog[1:n_events, :])
    data = getfield(trimmed_catalog, Symbol(quantity))

    stem!(ax, trimmed_catalog.t/(365*24*60*60), data; stem_kwargs...)

    return ax
end


function plotCatalogSSH(catalog_path::String, quantity::String; ax_kwargs, stem_kwargs)
    url = ENV["elja_url"]
    username = ENV["elja_user"]
    private_file = ENV["elja_private"]
    public_file = ENV["elja_pub"]


    sftp = SFTP(url, username, public_file, private_file)

    catalog  = load(SFTPClient.download(sftp, catalog_path))["data"]

    return plotCatalog(catalog, quantity; ax_kwargs, stem_kwargs)


end
function plotCatalogSSH(catalog_path::String, quantity::String, ax::Axis; stem_kwargs)
    url = ENV["elja_url"]
    username = ENV["elja_user"]
    private_file = ENV["elja_private"]
    public_file = ENV["elja_pub"]


    sftp = SFTP(url, username, public_file, private_file)

    catalog  = load(SFTPClient.download(sftp, catalog_path))["catalog"][1:upto, :]

    return plotCatalog(catalog, quantity, ax; stem_kwargs)


end
