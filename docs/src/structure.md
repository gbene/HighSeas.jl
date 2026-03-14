# Code structure and design


In this page the code structure and design are explained


I'm writing this mainly for myself so I don't get lost. I will describe how the package is structured at the moment.
I am trying to follow a composition based approach so that It is easier to customize the different solvers and problems but have a unified interface to solve/save/plot the results. 

## Base  structs
The following are structs that do not depend on anything:

### 1. AbstractMaterial
The structs that inherit from this type define the material of the experiment.  For example a simple material would be:

```julia
struct SimpleMaterial <: AbstractMaterial 
    cs::Float64
    rho::Float64
    nu::Float64
    a::Float64
    aRs::Float64
    b::Float64
    fr::Float64
    G::Float64
    eta::Float64
    Dc::Float64
end
```
### 2. AbstractGrid
The structs that inherit from this type define the base grid of the experiment.  For example a simple grid would be

```julia
struct Grid <: AbstractGrid
    cell_size::Float64
    W::Float64    # length of the entire domain
    L::Float64    # width of the entire domain

    n_elementsx::Int # number of elements on the x dimension
    n_elementsy::Int # number of elements on the y dimension

    X::StepRangeLen  # X range
    Y::StepRangeLen  # Y range

    x::AbstractArray{Float64} # x coordinates mesh
    y::AbstractArray{Float64} # y coordinates mesh

end
```


We can also make subtypes of AbstractGrid, for example a  `AbstractPowerGrid` that has specific needs (i.e. grid nodes **must** be a power of 2)

```julia
struct PowerGrid <: AbstractPowerGrid 
    
    cell_size::Float64
    W::Float64    # length of the entire domain
    L::Float64    # width of the entire domain


    domain_powerx::Int
    domain_powery::Int

    n_elementsx::Int # number of elements on the x dimension
    n_elementsy::Int # number of elements on the y dimension


    X::StepRangeLen  # X range
    Y::StepRangeLen  # Y range

    x::AbstractArray{Float64} # x coordinates mesh
    y::AbstractArray{Float64} # y coordinates mesh
end
```
### 3. AbstractState
The structs that inherit from this type define the state of the simulation (i.e. the variables that are being simulated).  For example a simple state would be:

```julia
struct State <: AbstractState
    
    dx::AbstractArray{Float64}
    V::AbstractArray{Float64}
    theta::AbstractArray{Float64}
    tau::AbstractArray{Float64}
end
```



## Constructed structs

These structs are constructed from base structs (i.e. use a struct to get info but don't have it as a field)


### 1. AbstractFault
The structs that inherit from this type define the geometry of the fault using the grid as reference. For example a rectangular fault would be

```julia
struct RectangleFault <: AbstractFault

    Wf::Float64 # half width of the fault
    Lf::Float64 # half length of the fault

    dLO::AbstractArray{Int8} # Mask indicating inside the fault
    dCR::AbstractArray{Int8} # Mask indicating outside the fault

    function RectangleFault(Wf::Float64, Lf::Float64, grid::AbstractGrid)

        x = grid.x
        y = grid.y
        

        dLO = @. Int8( (-Lf <= x <= Lf)*(-Wf <= y <= Wf));

        dCR = @. Int8(dLO == 0); 

        new(Wf, Lf, dLO, dCR)
    end

end
```

### 2. AbstractPatch
The structs that inherit from this type define the geometry of the slip patch (= rate weakening patch) on the fault using the grid as reference. For example a rectangular patch would be.

```julia
struct RectanglePatch <: AbstractPatch
    w::Float64  # half width of the rate weakening zone
    l::Float64  # half length of the rate weakening zone
    h::Float64  # Buffer between RW and RS

    dRW::AbstractArray{Int8} # Mask indicating inside the rate weakening zone
    dRS::AbstractArray{Int8} # Mask indicating inside the rate strengthening zone
    dTR::AbstractArray{Int8} # Mask indicating the transition between RS and RW

    function RectanglePatch(w::Float64, l::Float64, h::Float64, grid::AbstractGrid)
        x = grid.x
        y = grid.y

        dRW = @. Int8((abs(x) < l) * (abs(y) < w)); #inside rate/velocity weakening area
        dRS = @. Int8(dRW==0); #inside rate/velocity strengthening area
        dTR = @. Int8(((( l <= abs(x) <= l+h)*(abs(y)<=w+h))+(( w <= abs(y) <= w+h)*(abs(x)<=l+h)))>=1); #Transition zone between RS and RW


        new(w, l, h, dRW, dRS, dTR)

    end
    
end
```

### 3. AbstractNucleation
The structs that inherit from this type define a nucleation zone on the fault using the grid as reference. For example a rectangular nucleation would be

```julia
struct RectangleNucleation <: AbstractNucleation

    xi::Float64 # x center of the nucleation zone
    yi::Float64 # y center of the nucleation zone

    wi::Float64 # width of the nucleation zone
    li::Float64 # length of the nucleation zone

    dNU::AbstractArray{Int8} # Mask indicating inside the nucleation zone
    dFD::AbstractArray{Int8} # Mask indicationg everything outside the nucleation zone


    function RectangleNucleation(xi::Float64, yi::Float64, wi::Float64, li::Float64, grid::AbstractGrid)
        x = grid.x
        y = grid.y

        x_rlim = xi+li/2
        x_llim = xi-li/2

        y_ulim = yi+wi/2
        y_llim = yi-wi/2


        dNU = @. Int8((x_llim <= x < x_rlim)*(y_llim <= y < y_ulim))

        dFD = @. Int8(dNU==0)


        new(xi, yi, wi, li, dNU, dFD)

    end

end
```
## Depending structs

These structs have in the field other structs (either constructed and/or base). They are intended more as a means to group things together. 

### 1. AbstractDomain

The structs that inherit from this type define the entire domain that is being simulated. 

```julia
struct Domain <:AbstractDomain


    grid::AbstractGrid
    fault::AbstractFault
    patch::AbstractPatch
    nucleation::AbstractNucleation

    function Domain(grid::AbstractGrid, fault::AbstractFault, patch::AbstractPatch)

        nucleation = EmptyNucleation(grid)

        new(grid, fault, patch, nucleation)
    end

    function Domain(grid::AbstractGrid, fault::AbstractFault, patch::AbstractPatch, nucleation::AbstractNucleation)
        new(grid, fault, patch, nucleation)
    end

end

```

For example a pseudo code for a rectangular fault and patch build on a power grid would be 
```julia

grid = PowerGrid()

fault = RectangleFault()

patch = RectanglePatch()

domain = Domain(grid, fault, patch)
```

### 2. AbstractExperiment
The structs that inherit from this type define the experiment that we want to carry. This is where we do the setup of the experiment i.e. set initial values, a-b values, prestresses, check for the length scales etc etc. For example the BP4QD experiment is as follows:


<details>
  <summary>Spoiler warning</summary>

```julia
struct BP4QDExp <: AbstractBenchExperiment

    material::AbstractMaterial
    domain::AbstractDomain

    AL::Int # not sure if this goes here
    Vpl::Float64
    Vr::Float64
    Vi::Float64
    Vnu::Float64
    si::Float64

    lengthscales::NamedTuple

    template::AbstractArray{Float64} # not sure if this is necessary or not
    a::AbstractArray{Float64}
    b::AbstractArray{Float64}
    tau0::AbstractArray{Float64}

    state::AbstractState


    function BP4QDExp(input_dict::Dict, material::AbstractMaterial, domain::AbstractDomain)
        AL                = Int(input_dict["AL"])

        Vpl               = input_dict["Vpl"]
        Vi                = input_dict["Vi"]
        Vr                = input_dict["Vr"]
        Vnu               = input_dict["Vnu"]
        si                = input_dict["si0"]
        Dc                = input_dict["Dc"]


        grid = domain.grid
        patch = domain.patch
        nucleation = domain.nucleation

        aRSₛ = get_aRS(material)
        aₛ   = get_a(material)
        bₛ   = get_b(material)
        eta  = get_eta(material)
        fr   = get_fr(material)


        lengthscales = CheckLengthScales(material, domain, Dc, si)

        println(lengthscales)

        x = grid.x
        y = grid.y
        dRS = patch.dRS .== 1
        dTR = patch.dTR .== 1
        dNU = nucleation.dNU .== 1
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

        tau0 = @. si*a*asinh((Vi/(2*Vr))*exp((fr+bₛ*log(Vr/Vi))/a))+eta*Vi;
        @. tau0[dNU] = si*aₛ*asinh((Vnu/(2*Vr))*exp((fr+bₛ*log(Vr/Vi))/aₛ))+eta*Vnu;

        dx_init = copy(template)
        V_init  = Vi.+dx_init
        V_init[dNU] .= Vnu
        theta_init  = Dc/Vi.+dx_init

        tau_init = copy(tau0)

        state_init = State(dx_init, V_init, theta_init, tau_init)

        new(material, domain, AL, Vpl, Vi, Vr, Vnu, si, lengthscales, template, a, b, tau0, state_init)

    end

    function BP4QDExp(input_dict::Dict, material::AbstractMaterial, domain::AbstractDomain, state::AbstractState)

        Vpl               = input_dict["Vpl"]
        Vi                = input_dict["Vi"]
        Vr                = input_dict["Vr"]
        Vnu               = input_dict["Vnu"]
        si                = input_dict["si0"]


        grid = domain.grid
        patch = domain.patch
        nucleation = domain.nucleation

        aRSₛ = get_aRS(material)
        aₛ   = get_a(material)
        bₛ   = get_b(material)
        eta  = get_eta(material)
        fr   = get_fr(material)

        x = grid.x
        y = grid.y
        dRS = patch.dRS .== 1
        dTR = patch.dTR .== 1
        dNU = nucleation.dNU .== 1
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

        tau0 = @. si*a*asinh((Vi/(2*Vr))*exp((fr+bₛ*log(Vr/Vi))/a))+eta*Vi;
        @. tau0[dNU] = si*aₛ*asinh((Vnu/(2*Vr))*exp((fr+bₛ*log(Vr/Vi))/aₛ))+eta*Vnu;


        new(material, domain, Vpl, Vi, Vr, Vnu, si, template, a, b, tau0, state)

    end

end 
```

</details>


## Laws 
These type of structs are special because they are intended to be used as _functors_. This is because most of the laws require to pre allocate arrays for temporary or final results or for other things (such as FFT). Using the functor approach we can define a law as a struct and then use the object itself as a function without the need to pass the temporary matrices as input! This makes the code much more clean and general as each law will have references to the arrays/variables that are only necessary for their own calculations.

### AbstractLaw
The structs that inherit from this type define the laws that the experiment will obey. Each law will need to have a subtype defined so that in case we need a more granular control we can impose it. 

For example the AgeingLaw is a of type AbstractStateLaw that is a subtype ofAbstractLaw:

```julia
struct AgeingLaw <: AbstractStateLaw

    Dc::Float64
    theta:: AbstractArray{Float64}

    function AgeingLaw(experiment::AbstractExperiment)

        Dc = experiment.material.Dc
        theta = experiment.state.theta

        new(Dc, theta)
    end

end

# Define the functor. Only one per struct!!!
function (ageingLaw::AgeingLaw)(theta, V, dt)
    Dc = ageingLaw.Dc
    theta_res = ageingLaw.theta

    @. theta_res = theta * exp(-V * dt / Dc) + Dc / V * (1. - exp(-V * dt / Dc))

    return theta_res # this is optional but is nice to also have an output 
end
```

## What is missing?

- The final part is to put all these things together and define a type that is used to solve the problem. 
- Have a saver option
- We already have a stepper type that takes care of the time stepping, I did not describe it because the structure could change since it is tightly knit to the solver/algorithm part
- 


