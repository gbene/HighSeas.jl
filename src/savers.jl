struct EmptySaver <: AbstractSaver
end


"""
    StepSaver <: AbstractSaver


Save the state (dx, V, theta and tau) of the simulation at the given step.

### Fields

- `outpath::String` -- Where to save the output
- `every::Int` -- Save every nth step
- `data::AbstractState` -- State of the simulation to save
- `stepper::AbstractStepper` -- Stepper used to regulate the simulation


### Examples

- `StepSaver(every::Int, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)`

### Notes

- When using AbstractDetectors, `every` is ignored and the step is saved at the end of each detection
- When defining the saver, `outpath` is set from the `experiment`

"""
struct StepSaver{S<:AbstractState, ST<:AbstractStepper} <: AbstractSaver

    outpath::String
    every::Int
    data::S
    stepper::ST

    function StepSaver(every::Int, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

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

"""
    CatalogSaver <: AbstractSaver


Save the catalog of the simulation at the given step.

### Fields

- `outpath::String` -- Where to save the output
- `every::Int` -- Save every nth step
- `data::AbstractCatalog` -- Catalog of the simulation to save
- `stepper::AbstractStepper` -- Stepper used to regulate the simulation


### Examples

- `CatalogSaver(every::Int, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)`

### Notes

- When using AbstractDetectors, `every` is ignored and the step is saved at the end of each detection
- When defining the saver, `outpath` is set from the `experiment`

"""
struct CatalogSaver{C<:AbstractCatalog, ST<:AbstractStepper}  <: AbstractSaver
    outpath::String
    every::Int
    data::C
    stepper::ST

    function CatalogSaver(every::Int, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

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


"""
    SamplerSaver <: AbstractSaver


Save the samplers used during the simulation at the given step.

### Fields

- `outpath::String` -- Where to save the output
- `every::Int` -- Save every nth step
- `data::Vector{<:AbstractSampler}` -- Samplers to save
- `stepper::AbstractStepper` -- Stepper used to regulate the simulation


### Examples

- `SamplerSaver(every::Int, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)`

### Notes

- When using AbstractDetectors, `every` is ignored and the step is saved at the end of each detection
- When defining the saver, `outpath` is set from the `experiment`

"""
struct SamplerSaver{V<:Vector{<:AbstractSampler}, ST<:AbstractStepper} <: AbstractSaver
    outpath::String
    every::Int
    data::V
    stepper::ST

    function SamplerSaver(every::Int, samplers::Vector{<:AbstractSampler}, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

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


"""
    SnaptshotSaver <: AbstractSaver


Save AbstractPlotters during the simulation at the given step.

### Fields

- `outpath::String` -- Where to save the output
- `every::Int` -- Save every nth step
- `data::AbstractPlotter` -- State of the simulation to save
- `stepper::AbstractStepper` -- Stepper used to regulate the simulation


### Examples

- `SnaptshotSaver(every::Int, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)`

### Notes

- When using AbstractDetectors, `every` is ignored and the step is saved at the end of each detection
- When defining the saver, `outpath` is set from the `experiment`

"""
struct SnaptshotSaver{P<:AbstractPlotter, ST<:AbstractStepper} <: AbstractSaver

    outpath::String
    every::Int
    data::P
    stepper::ST

    function SnaptshotSaver(every::Int, plotter::AbstractPlotter, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

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
function simsave(saver::AbstractSaver, step::Int=1)

    if ~(typeof(saver) <: EmptySaver)

        if step == 1 || step % saver.every == 0
            filename, data, step = get_data(saver)

            @save filename data

        end
    end
end

function simsave(saver::SnaptshotSaver, step::Int=1)

    if step == 1 || step % saver.every == 0
        filename, data = get_data(saver)

        save(filename, data, px_per_unit=3)
    end
end


function simsave(savers::Vector{<:AbstractSaver})

    for saver in savers
        simsave(saver)
    end
end
