abstract type AbstractDetector end


function start_event_message(eventN, time, file)
    msg = "$(string(now())) event $eventN has started, time: $(time/(365*24*60*60))"
    println(msg)

    open(file,"a") do f
        write(f, "\n$msg")
    end

end
function end_event_message(mag, file)
    msg = "$(string(now())) Event has ended, mag: $mag"
    println(msg)

    open(file,"a") do f
        write(f, "\n $msg")
    end

end

function end_event_message(file)
    msg = "$(string(now())) Event has ended"
    println(msg)

    open(file,"a") do f
        write(f, "\n$msg")
    end

end



struct EmptyDetector <: AbstractDetector

end


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

    function CatalogDetector(minVThresh::Float64, maxVThresh::Float64, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)

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
            dx_start            = memcopy(dx_start)
            tau_start           = memcopy(tau_start)
            slip                = memcopy(slip)
            stressdrop          = memcopy(stressdrop)
            x                   = memcopy(x)
            y                   = memcopy(y)
            ruptured_nodes      = memcopy(ruptured_nodes)
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
