
mutable struct AdaptiveStepper{EL<:AbstractErrorLaw} <: AbstractAdaptiveStepper

    errorlaw::EL

    err::Float64
    erra::Float64

    tollo::Float64
    tolup::Float64
    mindt::Float64
    mintries::Int

    step::Int
    dt::Float64
    time::Float64
    breakcode::Bool


    function AdaptiveStepper(input_dict::Dict, errorlaw::AbstractErrorLaw, mintries::Int=2)


        frac = 1/2^(input_dict["fract"]-1)
        tollo = frac*input_dict["tollo"]
        tolup = frac*input_dict["tolup"]
        mindt = CalcMinDt(input_dict)
        mintries = mintries


        step = 0
        dt   = mindt
        time = 0.0
        breakcode=false

        err  = 1.0
        erra = 1.0

        new{typeof(errorlaw)}(errorlaw, err, erra, tollo, tolup, mindt, mintries, step, dt, time, breakcode)
    end

    function AdaptiveStepper(input_dict::Dict, errorlaw::AbstractErrorLaw, loadedstep::LoadedStep,  mintries::Int=2)


        frac = 1/2^(input_dict["fract"]-1)
        tollo = frac*input_dict["tollo"]
        tolup = frac*input_dict["tolup"]
        mindt = CalcMinDt(input_dict)
        mintries = mintries


        step = loadedstep.step
        dt   = mindt
        time = loadedstep.time
        breakcode=false

        err  = 1.0
        erra = 1.0

        new{typeof(errorlaw)}(errorlaw, err, erra, tollo, tolup, mindt, mintries, step, dt, time, breakcode)
    end

    function AdaptiveStepper{EL}(errorlaw::EL, err, erra, tollo, tolup, mindt, mintries, step, dt, time, breakcode) where EL
        new{EL}(errorlaw, err, erra, tollo, tolup, mindt, mintries, step, dt, time, breakcode)
    end
end


"""
Increase dt
"""
function StepUp(stepper::AbstractAdaptiveStepper)

    if stepper.err < stepper.tollo && stepper.erra < 10*stepper.tollo
        stepper.dt *= 1.1;
    end

    if stepper.dt < stepper.mindt
        stepper.dt = stepper.mindt;
    end

    return stepper.dt, stepper.breakcode

end

"""
Decrease dt
"""
function StepDown(stepper::AbstractAdaptiveStepper, ntries::Int, A, B)
        err = 1.0
        erra = 1.0
        errorlaw = stepper.errorlaw

        if ntries >= stepper.mintries
            err, erra = errorlaw(A, B);
        end


        if (ntries < stepper.mintries || ((err > stepper.tolup || erra > 10*stepper.tolup) && stepper.dt > stepper.mindt))
            if ntries > 1

                    # If accuracy metrics are violated then time-step is refined
                    stepper.dt /= 2.0;
            end
        else
            stepper.breakcode = true
        end

        stepper.err = err
        stepper.erra = erra

        return stepper.dt, stepper.breakcode
end

"""
Go to the next step == increase time and step
"""
function NextStep(stepper::AbstractStepper)
    stepper.breakcode = false
    stepper.time += stepper.dt
    stepper.step += 1
end

"""
Get current step info
"""
function GetStep(stepper::AbstractStepper)
    return (time=stepper.time, dt=stepper.dt, step=stepper.step)
end
