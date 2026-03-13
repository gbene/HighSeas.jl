
"""
    +(x::AbstractCatalog...)

Concatenate catalogs.

It is possible to concatenate and merge multiple catalogs by using the + operator


### Notes

This concatenation operator creates a new catalog without all zeros or all nans rows

### Examples
```julia

...

catalog = catalog_a+catalog_b+catalog_c

```
"""
function Base.:+(x::AbstractCatalog...)
    catalog = copy(x[1].catalog)
    for i in 2:length(x)
        c = x[i].catalog
        catalog = vcat(catalog, c)
    end

    mask = all.(!iszero, eachrow(catalog)) .* all.(!isnan, eachrow(catalog))


    return Catalog(catalog[mask, :])
end


"""
    State{AbstractArray{Float64}} <: AbstractState

State of the simulation at a given step at i.e. the values of dx, V, theta and tau

### Fields

- `dx::AbstractArray{Float64}` -- Slip at the step
- `V::AbstractArray{Float64}` -- Slip rate
- `theta::AbstractArray{Float64}` -- State at the step
- `tau::AbstractArray{Float64}` -- Shear stress at the step


### Notes

When creating an empty state the gpu_id can be specified to indicate on which GPU the data should reside.

### Examples

- `State(Ncols, Nrows; gpu_id=0)` -- Define a state of zeros for dx, V, theta and tau
"""
struct State{M<:AbstractArray{Float64}} <: AbstractState

    dx::M
    V::M
    theta::M
    tau::M

    function State(Ncols, Nrows; gpu_id::Int=0)

        dx     = zeros(Ncols, Nrows)
        V      = zeros(Ncols, Nrows)
        theta  = zeros(Ncols, Nrows)
        tau    = zeros(Ncols, Nrows)

        if typeof(get_backend()) <: AbstractGPUBackend
            dx      = memcopy(dx, gpu_id)
            V       = memcopy(V, gpu_id)
            theta   = memcopy(theta, gpu_id)
            tau     = memcopy(tau, gpu_id)
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


"""
    Catalog <: AbstractCatalog


Create a catalog.

This object contains all the event information calculated during the simulation

### Fields

- catalog::Matrix{Float64} -- Matrix representation of the catalog
- max_events::Int -- Maximum number of catalog events (i.e. number of rows)
- t::SubArray -- Time of rupter (in seconds)
- interevent_time::SubArray -- Interevent time between (in seconds)
- Moment::SubArray -- Moment released (in N⋅m)
- mag::SubArray -- Magnitude using Hanks and Kanamori 1979
- Area::SubArray -- Rupture area (m²)
- MeanSlip::SubArray -- Mean slip in the rupture area (m)
- MeanStress::SubArray -- Mean stress drop in the rupture area (Pa)
- hypo_x::SubArray -- x coordinate of the event (m)
- hypo_y::SubArray -- y coordinate of the event (m)
- n_events::Int -- Numbe of recorded events



### Bibliography

Hanks, Thomas C., and Hiroo Kanamori. ‘A Moment Magnitude Scale’. Journal of Geophysical Research: Solid Earth 84, no. B5 (1979): 2348–50. https://doi.org/10.1029/JB084iB05p02348.


"""
mutable struct Catalog{S<:SubArray} <: AbstractCatalog

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

        catalog = fill(NaN, max_events, 9)

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
