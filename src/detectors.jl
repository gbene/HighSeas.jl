abstract type AbstractDetector end



struct EmptyDetector <: AbstractDetector

end

function (emptyDetector::EmptyDetector)()

    if emptyDetector.eventStart == false
        emptyDetector.eventStart = true
    else
        emptyDetector.eventStart = false
    end
    return nothing
end

mutable struct SimpleDetector{S<:AbstractState, ST<:AbstractStepper} <: AbstractDetector

    eventStart::Bool
    eventN::Int
    maxVThresh::Float64
    state::S
    stepper::ST


    function SimpleDetector(maxVThresh::Float64, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)
        state = experiment.state
        stepper = algorithm.stepper
        new{typeof(state), typeof(stepper)}(false, 1, maxVThresh, state, stepper)
    end
end

function (simpleDetector::SimpleDetector)()

    if simpleDetector.eventStart == false
        simpleDetector.eventStart = true
    else
        simpleDetector.eventStart = false
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
    temp::M

    ruptured_nodes::B

    time_start::Float64
    last_event_time::Float64
    cell_area::Float64

    function CatalogDetector(minVThresh::Float64, maxVThresh::Float64, experiment::AbstractExperiment, algorithm::AbstractAlgorithm)
        Nx = experiment.domain.grid.n_elementsx
        Ny = experiment.domain.grid.n_elementsy

        dx_start=zeros(Nx, Ny)
        slip=zeros(Nx, Ny)
        tau_start=zeros(Nx, Ny)
        stressdrop=zeros(Nx, Ny)
        temp = zeros(Nx, Ny)
        ruptured_nodes=zeros(Int8, Nx, Ny)

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
            slip                = memcopy(slip)
            tau_start           = memcopy(tau_start)
            stressdrop          = memcopy(stressdrop)
            ruptured_nodes      = memcopy(ruptured_nodes)
        end

        new{typeof(state), typeof(stepper), typeof(catalog), typeof(dx_start), typeof(ruptured_nodes)}(false, 1, minVThresh, maxVThresh, state, stepper, material,
                                                                                     catalog, dx_start, slip, tau_start, stressdrop, x, y, temp,
                                                                                     ruptured_nodes, time_start, last_event_time, cell_area)
    end
end
function (catalogDetector::CatalogDetector)()

    state = catalogDetector.state
    stepper = catalogDetector.stepper
    material = catalogDetector.material
    catalog = catalogDetector.catalog
    eventN = catalogDetector.eventN

    temp = catalogDetector.temp

    if catalogDetector.eventStart == false
        catalogDetector.eventStart = true

        t = stepper.time

        copy!(catalogDetector.dx_start, state.dx)
        copy!(catalogDetector.tau_start, state.tau)

        catalogDetector.last_event_time = catalogDetector.time_start
        catalogDetector.time_start = t

        @.. thread=true catalogDetector.ruptured_nodes = state.V > catalogDetector.maxVThresh
        total_nodes = sum(catalogDetector.ruptured_nodes)

        catalog.t[eventN] = t

        @.. thread=true temp = catalogDetector.x * catalogDetector.ruptured_nodes
        catalog.hypo_x[eventN] = sum(temp)/total_nodes

        @.. thread=true temp = catalogDetector.y * catalogDetector.ruptured_nodes
        catalog.hypo_y[eventN] = sum(temp)/total_nodes



    else
        catalogDetector.eventStart = false

        trup = stepper.time-catalogDetector.time_start

        interevent_time = catalogDetector.time_start-catalogDetector.last_event_time

        @.. thread=true catalogDetector.slip = state.dx - catalogDetector.dx_start
        @.. thread=true catalogDetector.stressdrop = state.tau - catalogDetector.tau_start

        @.. thread=true catalogDetector.ruptured_nodes = catalogDetector.slip > material.Dc/2

        total_nodes = sum(catalogDetector.ruptured_nodes)

        Area = total_nodes*catalogDetector.cell_area

        @.. thread=true temp = catalogDetector.slip * catalogDetector.ruptured_nodes
        MeanSlip   = sum(temp)/total_nodes

        @.. thread=true temp = catalogDetector.stressdrop * catalogDetector.ruptured_nodes
        MeanStress = sum(temp)/total_nodes

        Moment = material.G*MeanSlip*Area
        mag = (log10(Moment)-9.05)/1.5

        catalog.interevent_time[eventN] = interevent_time
        catalog.Moment[eventN] = Moment
        catalog.mag[eventN] = mag
        catalog.Area[eventN] = Area
        catalog.MeanSlip[eventN] = MeanSlip
        catalog.MeanStress[eventN] = MeanStress

    end
    return nothing
end



function detect(detector::EmptyDetector)
    return nothing
end


function detect(detector::AbstractDetector)
    V = detector.state.V
    time = detector.stepper.time

    maxV = maximum(V)

    if maxV > detector.maxVThresh && detector.eventStart == false
            detector()
            println("event $(detector.eventN) has started, time: $(time/(365*24*60*60))")
    end

    if maxV <= detector.minVThresh && detector.eventStart == true
            detector()
            println("Event has ended")
            detector.eventN +=1
    end
    return nothing
end
