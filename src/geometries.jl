abstract type AbstractGrid end
abstract type AbstractPowerGrid <:AbstractGrid end

abstract type AbstractFault end
abstract type AbstractPatch end
abstract type AbstractNucleation end

abstract type AbstractDomain end


"""
Grid(input_dict::Dict)

# Description
Create a grid object.

# Attributes

+ cell_sizex::Float64
+ cell_sizey::Float64
+ cell_area::Float64
+ W::Float64    # length of the entire domain
+ L::Float64    # width of the entire domain
+ n_elementsx::Int # number of elements on the x dimension
+ n_elementsy::Int # number of elements on the y dimension
+ X::StepRangeLen  # X range
+ Y::StepRangeLen  # Y range
+ x::M # x coordinates mesh
+ y::M # y coordinates mesh

# Notes
"""
struct Grid{M<:AbstractArray{Float64}} <: AbstractGrid
    cell_sizex::Float64
    cell_sizey::Float64

    cell_area::Float64
    W::Float64    # length of the entire domain
    L::Float64    # width of the entire domain

    n_elementsx::Int # number of elements on the x dimension
    n_elementsy::Int # number of elements on the y dimension

    X::StepRangeLen  # X range
    Y::StepRangeLen  # Y range

    x::M # x coordinates mesh
    y::M # y coordinates mesh

    function Grid(input_dict::Dict)

        cell_sizex     = input_dict["cellsizex"]
        cell_sizey     = input_dict["cellsizey"]

        L             = input_dict["L"]
        W             = input_dict["W"]

        n_elementsx   = Int(L/cell_sizex)+1
        n_elementsy   = Int(W/cell_sizey)+1

        L_limits = L/2
        W_limits = W/2


        X = -L_limits:cell_sizex:L_limits
        Y = -W_limits:cell_sizey:W_limits

        x = X'.*ones(n_elementsy, n_elementsx) # The ys are the rows, xs are the columns
        y = Y.*ones(n_elementsy, n_elementsx)

        cell_area = cell_sizex*cell_sizey

        new{typeof(x)}(cell_sizex, cell_sizey, cell_area, L, W, n_elementsx, n_elementsy, X, Y, x, y)
    end

    function Grid{M}(cell_sizex, cell_sizey, cell_area, L, W, n_elementsx, n_elementsy, X, Y, x::M, y::M) where M
        new{M}(cell_sizex, cell_sizey, cell_area, L, W, n_elementsx, n_elementsy, X, Y, x, y)
    end


end


"""
PowerGrid(input_dict::Dict)

# Description
Create a grid object in which cells are square and the number of elements in the grid is a power of 2.
This means that a rounding operation is performed as follows round(log2(L/cell_size)), round(log2(W/cell_size))
to define the power of 2 for x and y.

# Attributes

+ cell_sizex::Float64
+ cell_sizey::Float64
+ cell_area::Float64
+ W::Float64                # length of the entire domain
+ L::Float64                # width of the entire domain
+ domain_powerx::Int        #Power of 2 along the x dimension
+ domain_powery::Int        #Power of 2 along the y dimension
+ n_elementsx::Int          # number of elements on the x dimension
+ n_elementsy::Int          # number of elements on the y dimension
+ X::StepRangeLen           # X range
+ Y::StepRangeLen           # Y range
+ x::M                      # x coordinates mesh
+ y::M                      # y coordinates mesh

# Notes
"""

struct PowerGrid{M<:AbstractArray{Float64}} <: AbstractPowerGrid

    cell_sizex::Float64
    cell_sizey::Float64

    cell_area::Float64

    W::Float64    # length of the entire domain
    L::Float64    # width of the entire domain


    domain_powerx::Int #Power of 2 along the x dimension
    domain_powery::Int #Power of 2 along the y dimension

    n_elementsx::Int # number of elements on the x dimension
    n_elementsy::Int # number of elements on the y dimension


    X::StepRangeLen  # X range
    Y::StepRangeLen  # Y range

    x::M # x coordinates mesh
    y::M # y coordinates mesh


    function PowerGrid(input_dict::Dict)

        cell_sizex     = input_dict["cellsizex"]
        cell_sizey     = input_dict["cellsizey"]

        @assert cell_sizex === cell_sizey "Cell is not square!"

        L             = input_dict["L"]
        W             = input_dict["W"]


        domain_powerx  = round(Int, log2(L/cell_sizex))
        domain_powery  = round(Int, log2(W/cell_sizey))

        n_elementsx   = 2^domain_powerx
        n_elementsy   = 2^domain_powery

        L_domain = (cell_sizex*n_elementsx)
        W_domain = (cell_sizey*n_elementsy)

        L_limits = L_domain/2
        W_limits = W_domain/2


        X = -L_limits:cell_sizex:L_limits-1
        Y = -W_limits:cell_sizey:W_limits-1

        x = X'.*ones(n_elementsy, n_elementsx) # The ys are the rows, xs are the columns
        y = Y.*ones(n_elementsy, n_elementsx)

        cell_area = cell_sizex*cell_sizey

        new{typeof(x)}(cell_sizex, cell_sizey, cell_area, L_domain, W_domain, domain_powerx, domain_powery, n_elementsx, n_elementsy, X, Y, x, y)
    end

    function PowerGrid{M}(cell_sizex, cell_sizey, cell_area, L_domain, W_domain, domain_powerx, domain_powery, n_elementsx, n_elementsy, X, Y, x::M, y::M) where M
        new{M}(cell_sizex, cell_sizey, cell_area, L_domain, W_domain, domain_powerx, domain_powery, n_elementsx, n_elementsy, X, Y, x, y)
    end

end

"""
RectangleFault(input_dict, grid)

Create a rectangular fault

# Attributes

+ Wf::Float64 # half width of the fault
+ Lf::Float64 # half length of the fault
+ dLO::AbstractArray{Int8} # Mask indicating the loading area
+ dCR::AbstractArray{Int8} # Mask indicating the creeping region

# Notes

"""
struct RectangleFault{M<:AbstractArray{Int8}} <: AbstractFault # Try if AbstractArray{Bool} works with FastBroadcasting

    Wf::Float64 # half width of the fault
    Lf::Float64 # half length of the fault

    dLO::M # Mask indicating inside the loading area
    dCR::M # Mask indicating the creeping region


    function RectangleFault(input_dict::Dict, grid::AbstractGrid)

        Wf            = input_dict["Wf"]
        Lf            = input_dict["Lf"]


        x = grid.x
        y = grid.y


        dLO = @. Int8( (-Lf <= x <= Lf)*(-Wf <= y <= Wf));
        dCR = @. Int8(dLO == 0);

        if typeof(get_backend()) <: AbstractGPUBackend
            dLO = memcopy(dLO)
            dCR = memcopy(dCR)
        end

        new{typeof(dLO)}(Wf, Lf, dLO, dCR)
    end

    function RectangleFault{M}(Wf, Lf, dLO::M, dCR::M) where M
        new{M}(Wf, Lf, dLO, dCR)
    end

end


"""
RectangleFault(input_dict, grid)

Create a rectangular slip patch

# Attributes

+ w::Float64  # half width of the rate weakening zone
+ l::Float64  # half length of the rate weakening zone
+ h::Float64  # Buffer between RW and RS
+ dRW::AbstractArray{Int8} # Mask indicating the rate weakening zone
+ dRS::AbstractArray{Int8} # Mask indicating the rate strengthening zone
+ dTR::AbstractArray{Int8} # Mask indicating the transition between RS and RW

# Notes

The buffer can be set to 0, this entails that there is no transition zone between RW and RS.
"""
struct RectanglePatch{M<:AbstractArray{Int8}} <: AbstractPatch
    w::Float64  # half width of the rate weakening zone
    l::Float64  # half length of the rate weakening zone
    h::Float64  # Buffer between RW and RS

    dRW::M # Mask indicating inside the rate weakening zone
    dRS::M # Mask indicating inside the rate strengthening zone
    dTR::M # Mask indicating the transition between RS and RW


    function RectanglePatch(input_dict::Dict, grid::AbstractGrid)
        w             = input_dict["w"]
        l             = input_dict["l"]
        h             = input_dict["h"]


        x = grid.x
        y = grid.y

        dRW = @. Int8((-l <= x <= l) * (-w <= y <= w)); #inside rate/velocity weakening area
        dRS = @. Int8(dRW==0); #inside rate/velocity strengthening area
        dTR = @. Int8(((( l <= abs(x) <= l+h)*(abs(y)<=w+h))+(( w <= abs(y) <= w+h)*(abs(x)<=l+h)))>=1); #Transition zone between RS and RW

        if typeof(get_backend()) <: AbstractGPUBackend
            dRW = memcopy(dRW)
            dRS = memcopy(dRS)
            dTR = memcopy(dTR)

        end

        new{typeof(dRW)}(w, l, h, dRW, dRS, dTR)

    end

    function RectanglePatch{M}(w, l, h, dRW::M, dRS::M, dTR::M) where M
        new{typeof(dRW)}(w, l, h, dRW, dRS, dTR)
    end

end


"""
CustomPatch(dRW::AbstractArray{Bool}, w::Float64, l::Float64, grid)

Create a custom slip patch object. A bool mask must be used as input together with the half width and length.

# Attributes

+ w::Float64  # half width of the rate weakening zone
+ l::Float64  # half length of the rate weakening zone
+ h::Float64  # Buffer between RW and RS
+ dRW::AbstractArray{Int8} # Mask indicating the rate weakening zone
+ dRS::AbstractArray{Int8} # Mask indicating the rate strengthening zone
+ dTR::AbstractArray{Int8} # Mask indicating the transition between RS and RW

# Notes

The buffer is by default NaN at the moment and no transition zone is supported.
"""
struct CustomPatch{M<:AbstractArray{Int8}} <: AbstractPatch
    w::Float64  # half width of the rate weakening zone
    l::Float64  # half length of the rate weakening zone
    h::Float64  # Buffer between RW and RS

    dRW::M # Mask indicating inside the rate weakening zone
    dRS::M # Mask indicating inside the rate strengthening zone
    dTR::M # Mask indicating the transition between RS and RW


    function CustomPatch(points::AbstractMatrix{Float64}, grid::AbstractGrid; w::Float64=NaN, l::Float64=NaN, h::Float64=NaN)


        x = grid.x
        y = grid.y

        g = [vec(x) vec(y)]
        shape = ClosedShape(points)





        mask = inpoly2(g, shape.points, shape.edges)[:,1]
        dRW = map(Int8, reshape(mask, size(x)))


        if isnan(h)
            dTR = zeros(Int8, grid.n_elementsy, grid.n_elementsx); #Transition zone between RS and RW
        else
            s = [1+2*h/shape.l, 1+2*h/shape.w]
            buffer = shape*s
            maskb = inpoly2(g, buffer.points, buffer.edges)[:,1]
            dTR = map(Int8, reshape(maskb, size(x)))
            dTR .-= dRW
        end

        if isnan(w)
            w = shape.w/2
        end
        if isnan(l)
            l = shape.l/2
        end

        dRS = @. Int8(dRW==0); #inside rate/velocity strengthening area



        if typeof(get_backend()) <: AbstractGPUBackend
            dRW = memcopy(dRW)
            dRS = memcopy(dRS)
            dTR = memcopy(dTR)

        end

        new{typeof(dRW)}(w, l, h, dRW, dRS, dTR)

    end

    function CustomPatch(points::AbstractMatrix{Float64}, grid::AbstractGrid, buffer_points::AbstractMatrix{Float64}; w::Float64=NaN, l::Float64=NaN)


        x = grid.x
        y = grid.y

        g = [vec(x) vec(y)]
        shape = ClosedShape(points)
        buffer = ClosedShape(buffer_points)


        mask = inpoly2(g, shape.points, shape.edges)[:,1]
        dRW = map(Int8, reshape(mask, size(x)))

        maskb = inpoly2(g, buffer.points, buffer.edges)[:,1]
        dTR = map(Int8, reshape(maskb, size(x)))
        dTR .-= dRW

        if isnan(w)
            w = shape.w/2
        end
        if isnan(l)
            l = shape.l/2
        end
        h = (buffer.l/2 - l)

        dRS = @. Int8(dRW==0); #inside rate/velocity strengthening area



        if typeof(get_backend()) <: AbstractGPUBackend
            dRW = memcopy(dRW)
            dRS = memcopy(dRS)
            dTR = memcopy(dTR)

        end

        new{typeof(dRW)}(w, l, h, dRW, dRS, dTR)

    end

    function CustomPatch{M}(w, l, h, dRW::M, dRS::M, dTR::M) where M
        new{M}(w, l, h, dRW, dRS, dTR)
    end

end

"""
EmptyNucleation(grid)


Create an empty nucleation object i.e. no nucleation patches in the domain.
"""
struct EmptyNucleation{M<:AbstractArray{Int8}} <: AbstractNucleation

    xi::Float64 # x center of the nucleation zone
    yi::Float64 # y center of the nucleation zone

    wi::Float64 # width of the nucleation zone
    li::Float64 # length of the nucleation zone

    dNU::M # Mask indicating inside the nucleation zone
    dFD::M # Mask indicationg everything outside the nucleation zone


    function EmptyNucleation(grid::AbstractGrid)
        xi = NaN
        yi = NaN

        wi = NaN
        li = NaN

        x = grid.x

        dNU = zeros(Int8, size(x))
        dFD = @. Int8(dNU==0)

        if typeof(get_backend()) <: AbstractGPUBackend
            dNU = memcopy(dNU)
            dFD = memcopy(dFD)

        end

        new{typeof(dNU)}(xi, yi, wi, li, dNU, dFD)

    end

    function EmptyNucleation{M}(xi, yi, wi, li, dNU::M, dFD::M) where M
        new{M}(xi, yi, wi, li, dNU, dFD)
    end


end

"""
RectangleNucleation(input_dict::Dict, grid)

Create a rectangular nucleation object i.e. a rectangular (or square) nucleation patch in the domain.

# Attributes


+ xi::Float64 # x center of the nucleation zone
+ yi::Float64 # y center of the nucleation zone
+ wi::Float64 # width of the nucleation zone
+ li::Float64 # length of the nucleation zone
+ dNU::AbstractArray{Int8} # Mask indicating inside the nucleation zone
+ dFD::AbstractArray{Int8} # Mask indicationg everything outside the nucleation zone

# Notes

"""
struct RectangleNucleation{M<:AbstractArray{Int8}} <: AbstractNucleation

    xi::Float64 # x center of the nucleation zone
    yi::Float64 # y center of the nucleation zone

    wi::Float64 # width of the nucleation zone
    li::Float64 # length of the nucleation zone

    dNU::M # Mask indicating inside the nucleation zone
    dFD::M # Mask indicationg everything outside the nucleation zone


    function RectangleNucleation(input_dict::Dict, grid::AbstractGrid)
        xi = input_dict["xi"]
        yi = input_dict["yi"]

        wi = input_dict["wi"]
        li = input_dict["li"]

        x = grid.x
        y = grid.y

        x_rlim = xi+li/2
        x_llim = xi-li/2

        y_ulim = yi+wi/2
        y_llim = yi-wi/2


        dNU = @. Int8((x_llim <= x < x_rlim)*(y_llim <= y < y_ulim)) # To be exact it should be <= but I am debugging

        dFD = @. Int8(dNU==0)

        if typeof(get_backend()) <: AbstractGPUBackend
            dNU = memcopy(dNU)
            dFD = memcopy(dFD)
        end

        new{typeof(dNU)}(xi, yi, wi, li, dNU, dFD)

    end

    function RectangleNucleation{M}(xi, yi, wi, li, dNU::M, dFD::M) where M
        new{M}(xi, yi, wi, li, dNU, dFD)
    end


end


"""
Domain(grid, fault, patch)
Domain(grid, fault, patch, nucleation)


Create the domain object.

# Attributes


+ xi::Float64 # x center of the nucleation zone
+ yi::Float64 # y center of the nucleation zone
+ wi::Float64 # width of the nucleation zone
+ li::Float64 # length of the nucleation zone
+ dNU::AbstractArray{Int8} # Mask indicating inside the nucleation zone
+ dFD::AbstractArray{Int8} # Mask indicationg everything outside the nucleation zone

# Notes

If no nucleation object is given then it is assumed that there is no nucleation patch.

"""
struct Domain{G<:AbstractGrid, F<:AbstractFault, P<:AbstractPatch, N<:AbstractNucleation} <: AbstractDomain


    grid::G
    fault::F
    patch::P
    nucleation::N

    function Domain(grid::G, fault::F, patch::P) where {G, F, P}

        nucleation = EmptyNucleation(grid)

        new{typeof(grid), typeof(fault), typeof(patch), typeof(nucleation)}(grid, fault, patch, nucleation)
    end

    function Domain(grid::G, fault::F, patch::P, nucleation::N) where {G, F, P, N}

        new{typeof(grid), typeof(fault), typeof(patch), typeof(nucleation)}(grid, fault, patch, nucleation)
    end

    function Domain{G, F, P, N}(grid::G, fault::F, patch::P, nucleation::N) where {G, F, P, N}
        new{G, F, P, N}(grid, fault, patch, nucleation)
    end

end
