using HighSeas
using HighSeas.TimeZones
using CUDA # Change this to AMDGPU, METAL and so on to use those backends. METAL does not support Float64 so for now it is not supported
using JLD2
using GLMakie
using BenchmarkTools

struct BP4QDExpGradient{F<:AbstractArray{Float64}, 
                        M<:HighSeas.AbstractMaterial, 
                        D<:HighSeas.AbstractDomain, 
                        S<:HighSeas.AbstractState, 
                        C<: HighSeas.AbstractCatalog} <: HighSeas.AbstractBenchExperiment

    material::M
    domain::D

    start_time::String
    outpath::String
    Vpl::Float64
    Vr::Float64
    Vi::Float64
    Vnu::Float64

    lengthscales::NamedTuple

    a::F
    b::F
    tau0::F
    si0::F

    state::S
    catalog::C


    function BP4QDExpGradient(input_dict::Dict, material::M, domain::D, gradient::F, n_events, output_dir::String; gpu_id::Int=0) where {M, D, F}

        start_time = string(now())
        println("Experiment start time: $start_time")
        outpath = HighSeas.make_outdir(start_time, output_dir)

        Vpl               = input_dict["Vpl"]
        Vi                = input_dict["Vi"]
        Vr                = input_dict["Vr"]
        Vnu               = input_dict["Vnu"]
        Dc                = input_dict["Dc"]


        grid = domain.grid
        patch = domain.patch
        nucleation = domain.nucleation

        aRSₛ = material.aRs
        aₛ   = material.a
        bₛ   = material.b
        eta  = material.eta
        fr   = material.fr

        si_avg = maximum(gradient)
        lengthscales = HighSeas.CheckLengthScales(material, domain, si_avg)

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

        si0 = gradient

        tau0 = @. si0*a*asinh((Vi/(2*Vr))*exp((fr+bₛ*log(Vr/Vi))/a))+eta*Vi;
        @. tau0[dNU] = si0[dNU]*a[dNU]*asinh((Vnu/(2*Vr))*exp((fr+bₛ*log(Vr/Vi))/a[dNU]))+eta*Vnu;

        dx_init = copy(template)
        V_init  = Vi.+dx_init
        V_init[dNU] .= Vnu
        theta_init  = Dc/Vi.+dx_init

        tau_init = copy(tau0)


        if typeof(get_backend()) <: HighSeas.AbstractGPUBackend
            a               = HighSeas.memcopy(a, gpu_id)
            b               = HighSeas.memcopy(b, gpu_id)
            tau0            = HighSeas.memcopy(tau0, gpu_id)
            si0             = HighSeas.memcopy(si0, gpu_id)
            dx_init         = HighSeas.memcopy(dx_init, gpu_id)
            V_init          = HighSeas.memcopy(V_init, gpu_id)
            theta_init      = HighSeas.memcopy(theta_init, gpu_id)
            tau_init        = HighSeas.memcopy(tau_init, gpu_id)
        end

        state_init = State(dx_init, V_init, theta_init, tau_init)
        catalog_init = Catalog(n_events)

        open("$outpath/simulation.log","w") do file
                write(file, "Experiment start time: $start_time \n")
                write(file, "================================================================\n")
                for k in sort!(collect(keys(input_dict)))
                    value = input_dict[k]
                    value_string = "$value"
                    if (k == "W") | (k=="L")
                        grid_value = getproperty(grid, Symbol(k))
                        if value != grid_value
                            value_string = "$grid_value #target was $value"
                        end
                    end
                    write(file, "$k: $value_string\n")
                end
                write(file, "================================================================\n")
                write(file, "$(string(lengthscales))\n")
            end


        new{typeof(a), typeof(material), typeof(domain), typeof(state_init), typeof(catalog_init)}(material, domain, start_time, outpath,
                                                                                                   Vpl, Vr, Vi, Vnu,
                                                                                                   lengthscales, a, b, tau0, si0,
                                                                                                   state_init, catalog_init)

    end


    function BP4QDExpGradient(input_dict::Dict, material::M, domain::D, gradient::F, n_events, output_dir::String, loadedstep::HighSeas.LoadedStep; gpu_id::Int=0) where {M, D, F}

        start_time = string(now())
        println("Experiment start time: $start_time")
        outpath = make_outdir(start_time, output_dir)
        state = loadedstep.state

        Vpl               = input_dict["Vpl"]
        Vi                = input_dict["Vi"]
        Vr                = input_dict["Vr"]
        Vnu               = input_dict["Vnu"]

        grid = domain.grid
        patch = domain.patch
        nucleation = domain.nucleation

        aRSₛ = material.aRs
        aₛ   = material.a
        bₛ   = material.b
        eta  = material.eta
        fr   = material.fr


        si_avg = mean(gradient)
        lengthscales = CheckLengthScales(material, domain, si_avg)

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

        si0 = gradient

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

    function BP4QDExpGradient{F, M, D, S, C}(material::M, domain::D, start_time, outpath,
                        Vpl, Vr, Vi, Vnu, lengthscales,
                        a::F, b::F, tau0::F, si0::F,
                        state_init::S, catalog_init::C) where {F, M, D, S, C}

        new{F, M, D, S, C}(material, domain, start_time,outpath,
                           Vpl, Vr, Vi, Vnu,
                           lengthscales, a, b, tau0, si0,
                           state_init, catalog_init)
    end
end


set_GPUbackend() # To use UnifiedMemory add "unified" as the method argument. Default is DeviceMemory. AMDGPU ony supports DeviceMemory
input_dict = readSheet("BP4input.txt")

# Define the domain

grid = PowerGrid(input_dict)
fault = RectangleFault(input_dict, grid)

patch = RectanglePatch(input_dict, grid)
nucleation = RectangleNucleation(input_dict, grid)

domain = Domain(grid, fault, patch, nucleation)

# Define the material

material = SimpleMaterial(input_dict)



# Define the experiment
si = input_dict["si0"]
gradient = @. si + 1e3 * abs(grid.y) # simple gradient along x. 

HighSeas.plotGradient(domain, gradient, 0.3, figdisplay=true)

experiment = BP4QDExpGradient(input_dict, material, domain, gradient, 10, "BP4QD_out")
# # experiment = BP4QDExp(input_dict, material, domain, 5, "BP4QD_out", saved_step)

# Define the Govenring equations

stresslaw = StressFFT(experiment)
statelaw = AgeingLaw(experiment)

explicitlaw = ExplicitRate(experiment)
linearlaw = LinearizedRate(experiment)
hybridlaw = HybridRate(0.00001, explicitlaw, linearlaw)

governing_equations = GoverningEquations(hybridlaw, stresslaw, statelaw)


# Define the error law and stepper

errorlaw = DoubleError(experiment)
stepper  = AdaptiveStepper(input_dict, errorlaw)
# stepper  = AdaptiveStepper(input_dict, errorlaw, saved_step)



# Define the algorithm used to solve

algorithm = CustomNewtonSolver(experiment, governing_equations, stepper)

# Define the detector (i.e. what to do when V reaches a certain value)

detector = CatalogDetector(1e-3, 1e-2, experiment, algorithm)

# Additional live plotting. This works only with CUDA.UnifiedMemory

# plotter = RSPlotter(experiment, algorithm, 10)

# Additional samplers

# samplers = Array{HighSeas.AbstractSampler, 1}(undef, 14)
# for sp in 1:14
#     samplers[sp] = PointSampler("../BP4_sample_points.txt", sp, 700000, experiment)
# end
# samplers[15] = SectionSampler("dx", 0.0, "y", 700000, experiment)
# samplers[16] = ContourSampler("V", 1e-3, experiment)
# samplers = [ContourSampler(:V, 1e-3, experiment)]

# Define savers

stepsaver = StepSaver(500, experiment, algorithm) # save simulation state at every 500 steps
catalogsaver = CatalogSaver(500, experiment, algorithm) # save catalog at every 500 steps

# # snapshotsaver = SnaptshotSaver("BP4QD_out", 100, plotter, experiment, algorithm) # save a figure at every 500 steps
# samplersaver = SamplerSaver(500, samplers, experiment, algorithm) # save the samplers at every 500 steps


# savers = [stepsaver, catalogsaver]
# or
# savers = [stepsaver, catalogsaver, snapshotsaver]
#or
savers = [stepsaver, catalogsaver]
#or ...




# Define the solver (i.e. use definite steps or time?)
tf = input_dict["tf"]*(365*24*60*60)

# solver = StepSolver(1, savers, detector)


solver = TimeSolver(tf, savers, detector)

# or add samplers
solver = TimeSolver(tf, savers, detector, samplers)


# Solve
HighSeas.solve(experiment, algorithm, solver)

# HighSeas.solve(experiment, algorithm, solver, plotter)
