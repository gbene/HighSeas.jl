

"""
    start_event_message(eventN, time, file)

Commodity function used by the detectors to print an start event message

### Input

- eventN -- Number of the event
- time -- time of the event (in s)
- file -- Output file of where to append the message
"""
function start_event_message(eventN, time, file)
    msg = "$(string(now())) event $eventN has started, time: $(time/(365*24*60*60))"
    println(msg)

    open(file,"a") do f
        write(f, "\n$msg")
    end

end


"""
    end_event_message(file)
    end_event_message(mag, file)

Commodity function used by the detectors to print an end event message, optionally with magnitude info

### Input

- mag -- Magnitude of the event
- file -- Output file of where to append the message
"""
function end_event_message(mag, file)
    msg = "$(string(now())) Event has ended, mag: $mag"
    println(msg)

    open(file,"a") do f
        write(f, "\n$msg")
    end

end

function end_event_message(file)
    msg = "$(string(now())) Event has ended"
    println(msg)

    open(file,"a") do f
        write(f, "\n$msg")
    end

end


"""
    EmptyDetector <: AbstractDetector

Detector that does nothing

"""
struct EmptyDetector <: AbstractDetector

end


"""
    SimpleDetector{AbstractState, AbstractStepper} <: AbstractDetector

Simple detector that just prints a message and logs the event in the simulation.log file

### Fields

- eventStart::Bool -- Flag used to define if the event has started or not
- eventN::Int -- Event number
- minVThresh::Float64 -- Minimum slip rate threshold to define the **END** of the event
- maxVThresh::Float64 -- Maximum slip rate threshold to define the **START** of the event
- state::AbstractState -- State of the simulation
- stepper::AbstractStepper -- Stepper used in the simulation
- log_file::String -- Path of where to log possible console outputs

### Notes

- When used, the detector saves the `AbstractSavers` of the simulation at the end of each event.
- The dector has an associated functor that can be used to detect the event



### Examples

```julia

...

detector = SimpleDetector(minVThresh, maxVThresh, experiment, algorithm) # define the detector

...

detector(savers) # Detect a possible event given the thresholds

```

"""
mutable struct SimpleDetector{S<:AbstractState, ST<:AbstractStepper} <: AbstractDetector

    eventStart::Bool
    eventN::Int
    minVThresh::Float64
    maxVThresh::Float64
    state::S
    stepper::ST
    log_file::String


    function SimpleDetector(minVThresh::Float64, maxVThresh::Float64, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)
        state = experiment.state
        stepper = algorithm.stepper
        log_file = "$(experiment.outpath)/simulation.log"
        new{typeof(state), typeof(stepper)}(false, 1, minVThresh, maxVThresh, state, stepper, log_file)
    end
end

function (simpleDetector::SimpleDetector)(savers)

    if simpleDetector.eventStart == false
        start_event_message(simpleDetector.eventN, simpleDetector.stepper.time, simpleDetector.log_file)
        simpleDetector.eventStart = true
    else
        simpleDetector.eventStart = false
        end_event_message(simpleDetector.log_file)
        simsave(savers)

    end
    return nothing
end


"""
    CatalogDetector{AbstractState, AbstractStepper, AbstractCatalog, AbstractArray{Float64}, AbstractArray{Int8}} <: AbstractDetector

Detector used to detect an event and calculate the associated values and fill a catalog

### Fields

- `eventStart::Bool` -- Flag used to define if the event has started or not
- `eventN::Int` -- Event number
- `minVThresh::Float64` -- Minimum slip rate threshold to define the **END** of the event
- `maxVThresh::Float64` -- Maximum slip rate threshold to define the **START** of the event
- `state::AbstractState` -- State of the simulation
- `stepper::AbstractStepper` -- Stepper used during simulation
- `material::AbstractMaterial` -- Material of the simulation
- `catalog::AbstractCatalog` -- Catalog used to saved the calculated values
- `dx_start::AbstractArray{Float64}` -- Values of slip at the start of the event. These are pre-allocated temporary matrices used to calculate the different quantities in the catalog
- `tau_start::AbstractArray{Float64}` -- Values of τ at the start of the event. These are pre-allocated temporary matrices used to calculate the different quantities in the catalog
- `slip::AbstractArray{Float64}` -- Values of slip at the end of the event. These are pre-allocated temporary matrices used to calculate the different quantities in the catalog
- `stressdrop::AbstractArray{Float64}` -- Stress drop at the start of the event. These are pre-allocated temporary matrices used to calculate the different quantities in the catalog
- `x::AbstractArray{Float64}` -- x coords of the grid. These are pre-allocated temporary matrices used to calculate the different quantities in the catalog
- `y::AbstractArray{Float64}` -- y coords of the grid. These are pre-allocated temporary matrices used to calculate the different quantities in the catalog
- `ruptured_nodes::B` -- Mask of which nodes are being ruptured
- `time_start::Float64` -- Start of the rupture
- `last_event_time::Float64` -- last rupture time (to calculate the intertime event)
- `cell_area::Float64` -- Area of the cells
- `log_file::String` -- Path of where to log possible console outputs



### Notes

- When used, the detector saves the `AbstractSavers` of the simulation at the end of each event.
- The dector has an associated functor that can be used to detect the event
- When using GPU, it is possible to define where the temporary matrices reside using `gpu_id`


### Examples

```julia

...

detector = CatalogDetector(minVThresh, maxVThresh, experiment, algorithm; gpu_id::Int=0) # define the detector

...

detector(savers) # Detect a possible event given the thresholds

```
"""
mutable struct CatalogDetector{S<:AbstractState, ST<:AbstractStepper, C<:AbstractCatalog, M<:AbstractArray{Float64}, B<:AbstractArray{Int8}} <: AbstractDetector

    eventStart::Bool
    eventN::Int
    minVThresh::Float64
    maxVThresh::Float64



    state::S
    stepper::ST
    material::AbstractMaterial
    catalog::C

    dx_start::M
    tau_start::M
    slip::M
    stressdrop::M
    x::M
    y::M
    ruptured_nodes::B

    time_start::Float64
    last_event_time::Float64
    cell_area::Float64
    log_file::String

    function CatalogDetector(minVThresh::Float64, maxVThresh::Float64, experiment::AbstractExperiment, algorithm::AbstractAlgorithm; gpu_id::Int=0)

        sz = size(experiment.state.dx)

        dx_start=zeros(sz)
        slip=zeros(sz)
        tau_start=zeros(sz)
        stressdrop=zeros(sz)
        ruptured_nodes=zeros(Int8, sz)

        time_start=0.0
        last_event_time=0.0

        state = experiment.state
        stepper = algorithm.stepper
        material = experiment.material
        catalog = experiment.catalog


        x = experiment.domain.grid.x
        y = experiment.domain.grid.y
        cell_area = experiment.domain.grid.cell_area

        if typeof(get_backend()) <: AbstractGPUBackend
            dx_start            = memcopy(dx_start, gpu_id)
            tau_start           = memcopy(tau_start, gpu_id)
            slip                = memcopy(slip, gpu_id)
            stressdrop          = memcopy(stressdrop, gpu_id)
            x                   = memcopy(x, gpu_id)
            y                   = memcopy(y, gpu_id)
            ruptured_nodes      = memcopy(ruptured_nodes, gpu_id)
        end

        log_file = "$(experiment.outpath)/simulation.log"

        new{typeof(state), typeof(stepper),
            typeof(catalog), typeof(dx_start),
            typeof(ruptured_nodes)}(false, 1, minVThresh, maxVThresh, state, stepper, material,
                                    catalog, dx_start, tau_start, slip, stressdrop, x, y,
                                    ruptured_nodes, time_start, last_event_time, cell_area, log_file)
    end
end

function (catalogDetector::CatalogDetector)(savers)

    state = catalogDetector.state
    stepper = catalogDetector.stepper
    material = catalogDetector.material
    catalog = catalogDetector.catalog
    eventN = catalogDetector.eventN
    ruptured_nodes = catalogDetector.ruptured_nodes

    t = stepper.time

    if catalogDetector.eventStart == false
        start_event_message(eventN, t, catalogDetector.log_file)

        catalogDetector.eventStart = true


        copy!(catalogDetector.dx_start, state.dx)
        copy!(catalogDetector.tau_start, state.tau)

        catalogDetector.last_event_time = catalogDetector.time_start
        catalogDetector.time_start = t

        @.. thread=true ruptured_nodes = state.V > catalogDetector.maxVThresh
        total_nodes = sum(ruptured_nodes)

        catalog.t[eventN] = t

        # @.. thread=true temp = catalogDetector.x * catalogDetector.ruptured_nodes
        catalog.hypo_x[eventN] = dot(catalogDetector.x, ruptured_nodes)/total_nodes

        # @.. thread=true temp = catalogDetector.y * catalogDetector.ruptured_nodes
        catalog.hypo_y[eventN] = dot(catalogDetector.y, ruptured_nodes)/total_nodes



    else
        catalogDetector.eventStart = false

        trup = t-catalogDetector.time_start

        interevent_time = catalogDetector.time_start-catalogDetector.last_event_time

        @.. thread=true catalogDetector.slip = state.dx - catalogDetector.dx_start
        @.. thread=true catalogDetector.stressdrop = state.tau - catalogDetector.tau_start

        @.. thread=true ruptured_nodes = catalogDetector.slip > material.Dc/2

        total_nodes = sum(ruptured_nodes)

        Area = total_nodes*catalogDetector.cell_area

        # @.. thread=true temp = catalogDetector.slip * catalogDetector.ruptured_nodes
        MeanSlip   = dot(catalogDetector.slip, ruptured_nodes)/total_nodes

        # @.. thread=true temp = catalogDetector.stressdrop * catalogDetector.ruptured_nodes
        MeanStress = dot(catalogDetector.stressdrop, ruptured_nodes)/total_nodes

        Moment = material.G*MeanSlip*Area
        mag = (log10(Moment)-9.05)/1.5

        catalog.interevent_time[eventN] = interevent_time
        catalog.Moment[eventN] = Moment
        catalog.mag[eventN] = mag
        catalog.Area[eventN] = Area
        catalog.MeanSlip[eventN] = MeanSlip
        catalog.MeanStress[eventN] = MeanStress
        catalog.n_events += 1
        end_event_message(mag, catalogDetector.log_file)
        simsave(savers)


    end
    return nothing
end



function detect(detector::EmptyDetector)
    return nothing
end


function detect(detector::AbstractDetector, savers)
    V = detector.state.V

    maxV = maximum(V)

    if maxV > detector.maxVThresh && detector.eventStart == false
            detector(savers)
    end

    if maxV <= detector.minVThresh && detector.eventStart == true
            detector(savers)
            detector.eventN +=1
    end
    return nothing
end
