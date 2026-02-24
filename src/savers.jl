abstract type AbstractSaver end

struct EmptySaver <: AbstractSaver
end

struct StepSaver{S<:AbstractState, ST<:AbstractStepper} <: AbstractSaver

    outpath::String
    every::Int
    data::S
    stepper::ST

    function StepSaver(outdir::String, every::Int, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

        # start_time = experiment.start_time
        # platform = get_backend().device

        # outpath = "$(outdir)/$(start_time)/$(platform)"

        # if ~isdir(outpath)
        #     mkpath(outpath)
        # end
        outpath = experiment.outpath

        data = experiment.state
        stepper = algorithm.stepper

        new{typeof(data), typeof(stepper)}(outpath, every, data, stepper)
    end
end


struct CatalogSaver{C<:AbstractCatalog, ST<:AbstractStepper}  <: AbstractSaver
    outpath::String
    every::Int
    data::C
    stepper::ST

    function CatalogSaver(outdir::String, every::Int, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

        # start_time = experiment.start_time
        # platform = get_backend().device

        # outpath = "$(outdir)/$(start_time)/$(platform)"

        # if ~isdir(outpath)
        #     mkpath(outpath)
        # end
        outpath = experiment.outpath

        data = experiment.catalog
        stepper = algorithm.stepper

        new{typeof(data), typeof(stepper)}(outpath, every, data, stepper)
    end
end

struct SamplerSaver{V<:Vector{<:AbstractSampler}, ST<:AbstractStepper} <: AbstractSaver
    outpath::String
    every::Int
    data::V
    stepper::ST

    function SamplerSaver(outdir::String, every::Int, samplers::Vector{<:AbstractSampler}, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

        # start_time = experiment.start_time
        # platform = get_backend().device

        # outpath = "$(outdir)/$(start_time)/$(platform)"

        # if ~isdir(outpath)
        #     mkpath(outpath)
        # end
        outpath = experiment.outpath

        stepper = algorithm.stepper
        data = samplers

        new{typeof(data), typeof(stepper)}(outpath, every, data, stepper)

    end
end

struct SnaptshotSaver{P<:AbstractPlotter, ST<:AbstractStepper} <: AbstractSaver

    outpath::String
    every::Int
    data::P
    stepper::ST

    function SnaptshotSaver(outdir::String, every::Int, plotter::AbstractPlotter, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

        # start_time = experiment.start_time
        # platform = get_backend().device

        # outpath = "$(outdir)/$(start_time)/$(platform)"

        # if ~isdir(outpath)
        #     mkpath(outpath)
        # end
        outpath = "$(experiment.outpath)/saved_figures"

        # outpath = "$(outdir)/$(start_time)/$(device)/saved_figures"
        if ~isdir(outpath)
            mkpath(outpath)
        end

        stepper = algorithm.stepper

        # if get_backend().platform == "CPU"
        #     CairoMakie.activate!()
        # end

        data = plotter
        new{typeof(data), typeof(stepper)}(outpath, every, data, stepper)
    end
end

"""
Internal function used to get the necessary data of an AbstractSaver to save. Add type specific
functions for AbstractSavers that require specific treatment (e.g. the SnaptshotSaver)
"""
function get_data(saver::AbstractSaver)
    data = saver.data
    step = GetStep(saver.stepper)
    name = string(nameof(typeof(saver)))
    filename = "$(saver.outpath)/saved_$name.jld2"

    return filename, data, step
end

function get_data(saver::StepSaver)
    data = (saver.data, saver.stepper.step, saver.stepper.time)
    step = GetStep(saver.stepper)
    name = string(nameof(typeof(saver)))
    filename = "$(saver.outpath)/saved_$name.jld2"

    return filename, data, step
end


function get_data(saver::SnaptshotSaver)

    plotter = saver.data
    step = GetStep(saver.stepper)
    UpdatePlot(plotter)
    fig = plotter.fig

    if get_backend().platform == "CPU"
        filename = "$(saver.outpath)/$(step.step).svg"
    else
        filename = "$(saver.outpath)/$(step.step).png"
    end

    return filename, fig
end


"""
Common interface function to save any AbstractSaver. Add type specific functions
for AbstractSavers that require specific treatment (e.g. the SnaptshotSaver)
"""
function simsave(saver::AbstractSaver)

    if ~(typeof(saver) <: EmptySaver)
        step = saver.stepper.step


        if step % saver.every == 0 || step == 1
            filename, data, step = get_data(saver)

            @save filename data

        end
    end
end

function simsave(saver::SnaptshotSaver)

    step = saver.stepper.step

    if step % saver.every == 0 || step == 1
        filename, data = get_data(saver)

        save(filename, data, px_per_unit=3)
    end
end


function simsave(savers::Vector{<:AbstractSaver})

    for saver in savers
        simsave(saver)
    end
end
