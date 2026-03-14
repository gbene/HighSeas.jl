


function plot_experiment(experiment::AbstractExperiment, field::Symbol, inspector=false)
    fig = Figure(inspectable=inspector)
    ax = Axis(fig[1,1], aspect = DataAspect(), title = String(field))

    grid = experiment.domain.grid

    data = getfield(experiment, field)
    hm = heatmap!(ax, grid.X, grid.Y, data')
    Colorbar(fig[1,2], hm)
    DataInspector(fig)

    return fig, ax

end
function plot_experiment(experiment::AbstractExperiment, field::Symbol,  scale::Function, inspector=false)
    fig = Figure(inspectable=inspector)
    ax = Axis(fig[1,1], aspect = DataAspect(), title = String(field))

    grid = experiment.domain.grid

    data = getfield(experiment, field)
    hm = heatmap!(ax, grid.X, grid.Y, scale.(data'))
    Colorbar(fig[1,2], hm)
    DataInspector(fig)

    return fig, ax

end


"""
    BP4QDExp{AbstractArray{Float64}, AbstractMaterial, AbstractDomain, AbstractState, AbstractCatalog}

Create the experiment object for the BP4QD benchmark. N_events sets the size of the Catalog.

### Fields

+ material::AbstractMaterial -- Material of the experiment
+ domain::AbstractDomain -- Domain of the experiment
+ start_time::String -- Start time of the experiment
+ outpath::String -- Where to save the output of the simulations
+ Vpl::Float64 -- Velocity of the plate (m/s)
+ Vr::Float64 -- Reference velocity (m/s)
+ Vi::Float64 -- Initial velocity (m/s)
+ Vnu::Float64 -- Nucleation patch velocity (m/s)
+ lengthscales::NamedTuple -- Calculated lengthscales (Lb, Linf, L/Linf)
+ a::AbstractArray{Float64} -- Matrix of a values for rate and state
+ b::AbstractArray{Float64} -- Matrix of b values for rate and state
+ tau0::AbstractArray{Float64} -- Matrix of initial tau values for rate and state (Pa)
+ si0::AbstractArray{Float64} -- Matrix of initial si for rate and state (Pa)
+ state::AbstractState -- State of the simulation (this is not θ of the rate and state!!)
+ catalog::AbstractCatalog -- Catalog of the simulation


### Notes

- When defining the experiment from scratch (i.e. not making a copy) the simulation state (`::AbstractState`) and catalog (`::AbstractCatalog`) are automatically created
- When using GPUs, it is possible to set the `gpu_id` that defines where the temporary matrices reside
- If an experiment must be resumed by a backup, it is possible to pass a `LoadedStep` object

When the experiment is created, the following file structure is generated

```bash
.
└── Timestamp/
    └── Backend/
        └── simulation.log
```

Where the simulation.log will have the input sheet + eventual appended messages by other components (i.e. the `AbstractDetectors`).

### Examples

- `BP4QDExp(input_dict::Dict, material::AbstractMaterial, domain::AbstractDomain, n_events::Int; gpu_id=0)` -- Create a fresh experiment
- `BP4QDExp(input_dict, material, domain, n_events, output_dir, loadedstep::LoadedStep; gpu_id=0)` -- Resume an experiment from a given step

"""
struct BP4QDExp{F<:AbstractArray{Float64}, M<:AbstractMaterial, D<:AbstractDomain, S<:AbstractState, C<: AbstractCatalog} <: AbstractBenchExperiment

    material::M
    domain::D

    start_time::String
    outpath::String
    Vpl::Float64
    Vr::Float64
    Vi::Float64
    Vnu::Float64
    # si::Float64

    lengthscales::NamedTuple

    # template::M # not sure if this is necessary or not
    a::F
    b::F
    tau0::F
    si0::F


    # These are the values that are then updated at every step
    # I have put it here but I am still not super sure if it is the right place.
    # I was thinking on the lines of "imagine that we have an experimental setup that can measure these quantities"
    # i.e. here in AbstractExperiments

    # As of now I have added two init possibilities. The first that calculates the initial values
    # the second that gets the initial values. This will help I think in the future for resuming the simulation at a give step
    # Maybe intead of passing the single arrays we could pass a struct? Need to think about it

    # Update: Seems to be nice and clean. We can reference the single variables in the state struct as soft copies (pointing to the same memory)
    # Thus we can use them no problem in the different laws and it updates the values in the state struct.

    # dx::AbstractArray{Float64}
    # V::AbstractArray{Float64}
    # theta::AbstractArray{Float64}
    # tau::AbstractArray{Float64}
    # dt::Float64

    state::S
    catalog::C


    function BP4QDExp(input_dict::Dict, material::M, domain::D, n_events, output_dir::String; gpu_id::Int=0) where {M, D}

        start_time = string(now())
        println("Experiment start time: $start_time")
        outpath = make_outdir(start_time, output_dir)

        Vpl               = input_dict["Vpl"]
        Vi                = input_dict["Vi"]
        Vr                = input_dict["Vr"]
        Vnu               = input_dict["Vnu"]
        si                = input_dict["si0"]
        Dc                = input_dict["Dc"]


        grid = domain.grid
        patch = domain.patch
        nucleation = domain.nucleation

        aRSₛ = material.aRs
        aₛ   = material.a
        bₛ   = material.b
        eta  = material.eta
        fr   = material.fr


        lengthscales = CheckLengthScales(material, domain, si)

        println(lengthscales)

        x = grid.x
        y = grid.y
        dRS = Matrix(patch.dRS .== 1) # i really hate this, maybe we should have a function that before launching the solve moves everything to GPU?
        dTR = Matrix(patch.dTR .== 1)
        dNU = Matrix(nucleation.dNU .== 1)
        L = patch.l
        H = patch.w
        h = patch.h

        x_buffer = x[dTR]
        y_buffer = y[dTR]


        r = maximum([(abs.(x_buffer).- L);; (abs.(y_buffer).- H)],dims=2)./h;



        template = zeros(size(x))

        a_buffer = @. aₛ + r*(aRSₛ-aₛ);
        a = aₛ .+ template;
        a[dRS] .= aRSₛ;
        a[dTR] .= a_buffer;
        b = bₛ .+ template;

        si0 = si .+template

        tau0 = @. si0*a*asinh((Vi/(2*Vr))*exp((fr+bₛ*log(Vr/Vi))/a))+eta*Vi;
        @. tau0[dNU] = si0[dNU]*a[dNU]*asinh((Vnu/(2*Vr))*exp((fr+bₛ*log(Vr/Vi))/a[dNU]))+eta*Vnu;

        dx_init = copy(template)
        V_init  = Vi.+dx_init
        V_init[dNU] .= Vnu
        theta_init  = Dc/Vi.+dx_init

        tau_init = copy(tau0)


        if typeof(get_backend()) <: AbstractGPUBackend
            a               = memcopy(a, gpu_id)
            b               = memcopy(b, gpu_id)
            tau0            = memcopy(tau0, gpu_id)
            si0             = memcopy(si0, gpu_id)
            dx_init         = memcopy(dx_init, gpu_id)
            V_init          = memcopy(V_init, gpu_id)
            theta_init      = memcopy(theta_init, gpu_id)
            tau_init        = memcopy(tau_init, gpu_id)
        end

        state_init = State(dx_init, V_init, theta_init, tau_init)
        catalog_init = Catalog(n_events)

        open("$outpath/simulation.log","w") do file
                write(file, "Experiment start time: $start_time \n")
                write(file, "================================================================\n")
                for k in sort!(collect(keys(input_dict)))
                        write(file, "$k: $(input_dict[k])\n")
                end
                write(file, "================================================================\n")
                write(file, "$(string(lengthscales))\n")
            end


        new{typeof(a), typeof(material), typeof(domain), typeof(state_init), typeof(catalog_init)}(material, domain, start_time, outpath,
                                                                                                   Vpl, Vr, Vi, Vnu,
                                                                                                   lengthscales, a, b, tau0, si0,
                                                                                                   state_init, catalog_init)

    end


    function BP4QDExp(input_dict::Dict, material::M, domain::D, n_events, output_dir::String, loadedstep::LoadedStep; gpu_id::Int=0) where {M, D}

        start_time = string(now())
        println("Experiment start time: $start_time")
        outpath = make_outdir(start_time, output_dir)
        state = loadedstep.state

        Vpl               = input_dict["Vpl"]
        Vi                = input_dict["Vi"]
        Vr                = input_dict["Vr"]
        Vnu               = input_dict["Vnu"]
        si                = input_dict["si0"]
        Dc                = input_dict["Dc"]


        grid = domain.grid
        patch = domain.patch
        nucleation = domain.nucleation

        aRSₛ = material.aRs
        aₛ   = material.a
        bₛ   = material.b
        eta  = material.eta
        fr   = material.fr


        lengthscales = CheckLengthScales(material, domain, si)

        println(lengthscales)

        x = grid.x
        y = grid.y
        dRS = Matrix(patch.dRS .== 1) # i really hate this, maybe we should have a function that before launching the solve moves everything to GPU?
        dTR = Matrix(patch.dTR .== 1)
        dNU = Matrix(nucleation.dNU .== 1)
        L = patch.l
        H = patch.w
        h = patch.h

        x_buffer = x[dTR]
        y_buffer = y[dTR]


        r = maximum([(abs.(x_buffer).- L);; (abs.(y_buffer).- H)],dims=2)./h;



        template = zeros(size(x))

        a_buffer = @. aₛ + r*(aRSₛ-aₛ);
        a = aₛ .+ template;
        a[dRS] .= aRSₛ;
        a[dTR] .= a_buffer;
        b = bₛ .+ template;

        si0 = si .+template

        tau0 = @. si0*a*asinh((Vi/(2*Vr))*exp((fr+bₛ*log(Vr/Vi))/a))+eta*Vi;
        @. tau0[dNU] = si0[dNU]*a[dNU]*asinh((Vnu/(2*Vr))*exp((fr+bₛ*log(Vr/Vi))/a[dNU]))+eta*Vnu;

        dx_init = state.dx
        V_init  = state.V
        theta_init  = state.theta
        tau_init = state.tau


        if typeof(get_backend()) <: AbstractGPUBackend
            a               = memcopy(a, gpu_id)
            b               = memcopy(b, gpu_id)
            tau0            = memcopy(tau0, gpu_id)
            si0             = memcopy(si0, gpu_id)
        end

        state_init = State(dx_init, V_init, theta_init, tau_init)

        catalog_init = Catalog(n_events)

        open("$outpath/simulation.log","w") do file
            write(file, "Experiment start time: $start_time \n")
            write(file, "================================================================\n")
            for k in sort!(collect(keys(input_dict)))
                    write(file, "$k: $(input_dict[k])\n")
            end
            write(file, "================================================================\n")
            write(file, "$(string(lengthscales))\n")
        end


        new{typeof(a), typeof(material), typeof(domain), typeof(state_init), typeof(catalog_init)}(material, domain, start_time, outpath,
                                                                                                   Vpl, Vr, Vi, Vnu,
                                                                                                   lengthscales, a, b, tau0, si0,
                                                                                                   state_init, catalog_init)

    end

    function BP4QDExp{F, M, D, S, C}(material::M, domain::D, start_time, outpath,
                        Vpl, Vr, Vi, Vnu, lengthscales,
                        a::F, b::F, tau0::F, si0::F,
                        state_init::S, catalog_init::C) where {F, M, D, S, C}

        new{F, M, D, S, C}(material, domain, start_time,outpath,
                           Vpl, Vr, Vi, Vnu,
                           lengthscales, a, b, tau0, si0,
                           state_init, catalog_init)
    end
end
