


struct CustomNewtonSolver{E<:AbstractGoverningEquations, S<:AbstractAdaptiveStepper, I<:AbstractArray{Int8}, M<:AbstractArray{Float64}} <: AbstractNewton

    # stresslaw::AbstractStressLaw   # Law to calculate the stress
    # statelaw::AbstractStateLaw     # Law to calculate the state
    # ratelaw::AbstractRateLaw       # Law to calculate the rate
    ntries::Int
    equations::E
    stepper::S         # Stepper used to decide the step in the simulation


    dLO::I
    dCR::I
    Vpl::Float64

    si0::M
    tau0::M
    a::M
    b::M



    # Preallocated arrays to store temporary results

    dxp::M         # Slip at the step
    Vp::M          # Rate at the step
    thetap::M      # State at the step

    dxg::M        # Guess value of dx
    Vg::M         # Guess value of V
    Vmg::M        # Average between Vp and Vg
    Vm::M         # Average between Vp and V

    #thetag is not needed as preallocation, we can use the temporary value of the state!


    function CustomNewtonSolver{E, S, I, M}(ntries, equations::E, stepper::S, dLO::I, dCR::I,
                                     Vpl, si0::M, tau0::M, a::M, b::M,
                                     dxp::M, Vp::M, thetap::M,
                                     dxg::M, Vg::M, Vmg::M, Vm::M) where {E, S, I, M}


        new{E, S, I, M}(ntries, equations, stepper, dLO, dCR,
                        Vpl, si0, tau0, a, b, dxp,
                        Vp, thetap, dxg, Vg, Vmg, Vm)
    end

    function CustomNewtonSolver(experiment::AbstractExperiment, equations::AbstractGoverningEquations, stepper::AbstractStepper, ntries=10; gpu_id::Int=0)

        domainsize = size(experiment.domain.grid.x)

        dLO         = experiment.domain.fault.dLO # = is never a deep copy, it is always a reference to the original object
        dCR         = experiment.domain.fault.dCR
        Vpl         = experiment.Vpl


        si0         = experiment.si0
        tau0        = experiment.tau0
        a           = experiment.a
        b           = experiment.b



        dxp     = zeros(domainsize)
        Vp      = zeros(domainsize)
        thetap  = zeros(domainsize)

        dxg     = zeros(domainsize)
        Vg      = zeros(domainsize)
        Vmg     = zeros(domainsize)
        Vm      = zeros(domainsize)

        if typeof(get_backend()) <: AbstractGPUBackend
            dxp     = memcopy(dxp, gpu_id)
            Vp      = memcopy(Vp, gpu_id)
            thetap  = memcopy(thetap, gpu_id)
            dxg     = memcopy(dxg, gpu_id)
            Vg      = memcopy(Vg, gpu_id)
            Vmg     = memcopy(Vmg, gpu_id)
            Vm      = memcopy(Vm, gpu_id)
        end



        new{typeof(equations), typeof(stepper), typeof(dLO), typeof(si0)}(ntries, equations, stepper, dLO, dCR, Vpl,
                                                                          si0, tau0, a, b, dxp, Vp, thetap,
                                                                          dxg, Vg, Vmg, Vm)
    end

end



function (customNewtonSolver::CustomNewtonSolver)(dx, V, theta)

    dLO         = customNewtonSolver.dLO # this needs to be added to solver because it allocates and it takes longer to access
    dCR         = customNewtonSolver.dCR
    Vpl         = customNewtonSolver.Vpl

    si0         = customNewtonSolver.si0
    tau0        = customNewtonSolver.tau0
    a           = customNewtonSolver.a
    b           = customNewtonSolver.b


    dxp         = customNewtonSolver.dxp
    Vp          = customNewtonSolver.Vp
    thetap      = customNewtonSolver.thetap

    dxg         = customNewtonSolver.dxg
    Vg          = customNewtonSolver.Vg
    Vmg         = customNewtonSolver.Vmg

    Vm          = customNewtonSolver.Vm


    stresslaw   = customNewtonSolver.equations.stresslaw
    statelaw    = customNewtonSolver.equations.statelaw
    ratelaw     = customNewtonSolver.equations.ratelaw
    stepper     = customNewtonSolver.stepper

    n_tries     = customNewtonSolver.ntries


    copy!(dxp, dx)
    copy!(Vp, V)
    copy!(thetap, theta)

    dt, exitcode = StepUp(stepper)


    for trycount in 1:n_tries
        # If this is the first attempt
        if trycount == 1

            copy!(Vg, V);
            copy!(Vmg, V);
            @.. thread=true dxg = dx + dt*V;

        elseif trycount == 2
            # First iteration is treated as a special case where time-step
            # is not decreased, in attempt to obtain a better guess.

            copy!(dxg, dx);
            copy!(Vg, V);
            copy!(Vmg, Vm);

        else

            @.. thread=true dxg = 0.5*(dxp + dx);
            @.. thread=true Vg  = 0.5*(Vp + V);
            @.. thread=true Vmg = 0.5*(Vp + Vg);

        end

        tau = stresslaw(dxg, tau0)


        thetag = statelaw(thetap, Vmg, dt)


        V = ratelaw(Vg, tau, thetag, si0, a, b)


        @.. thread=true V = V*dLO + Vpl*dCR
        @.. thread=true Vm = 0.5*(Vp + V);
        @.. thread=true dx = dxp + dt*Vm;

        dt, exitcode = StepDown(stepper, trycount, V, Vg)


        if exitcode
            break
        end

    end

    theta = statelaw(thetap, Vm, dt)
    NextStep(stepper)
    return dx, V, theta


end
