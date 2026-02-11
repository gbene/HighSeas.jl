abstract type AbstractState end
abstract type AbstractCatalog end


struct State{M<:AbstractArray{Float64}} <: AbstractState

    dx::M
    V::M
    theta::M
    tau::M

    function State(Ncols, Nrows)

        dx     = zeros(Ncols, Nrows)
        V      = zeros(Ncols, Nrows)
        theta  = zeros(Ncols, Nrows)
        tau    = zeros(Ncols, Nrows)

        if typeof(get_backend()) <: AbstractGPUBackend
            dx      = memcopy(dx)
            V       = memcopy(V)
            theta   = memcopy(theta)
            tau     = memcopy(tau)
        end

        new{typeof(dx)}(dx, V, theta, tau)
    end
    function State(dx::M, V::M, theta::M, tau::M) where M
        new{typeof(dx)}(dx, V, theta, tau)
    end
    function State{M}(dx::M, V::M, theta::M, tau::M) where M
        new{M}(dx, V, theta, tau)
    end
end

struct Catalog{S<:SubArray} <: AbstractCatalog

    catalog::Matrix{Float64}
    max_events::Int

    t::S
    interevent_time::S
    Moment::S
    mag::S
    Area::S
    MeanSlip::S
    MeanStress::S
    hypo_x::S
    hypo_y::S
    n_events::Int



    function Catalog(max_events::Int)

        catalog = zeros(max_events, 9)

        t                   = @view catalog[:,1]
        interevent_time     = @view catalog[:,2]
        Moment              = @view catalog[:,3]
        mag                 = @view catalog[:,4]
        Area                = @view catalog[:,5]
        MeanSlip            = @view catalog[:,6]
        MeanStress          = @view catalog[:,7]
        hypo_x              = @view catalog[:,8]
        hypo_y              = @view catalog[:,9]
        n_events = 0

        new{typeof(t)}(catalog, max_events, t, interevent_time, Moment, mag, Area, MeanSlip, MeanStress, hypo_x, hypo_y, n_events)

    end

    function Catalog(catalog::Matrix{Float64})

        t                   = @view catalog[:,1]
        interevent_time     = @view catalog[:,2]
        Moment              = @view catalog[:,3]
        mag                 = @view catalog[:,4]
        Area                = @view catalog[:,5]
        MeanSlip            = @view catalog[:,6]
        MeanStress          = @view catalog[:,7]
        hypo_x              = @view catalog[:,8]
        hypo_y              = @view catalog[:,9]

        max_events = n_events = size(catalog)[1]

        new{typeof(t)}(catalog, max_events, t, interevent_time, Moment, mag, Area, MeanSlip, MeanStress, hypo_x, hypo_y, n_events)
    end
    function Catalog{S}(catalog, max_events, t::S, interevent_time::S,
                    Moment::S, mag::S, Area::S, MeanSlip::S,
                    MeanStress::S, hypo_x::S, hypo_y::S, n_events) where S

        new{S}(catalog, max_events, t, interevent_time, Moment, mag, Area, MeanSlip, MeanStress, hypo_x, hypo_y, n_events)
    end
end
