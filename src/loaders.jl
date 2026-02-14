abstract type AbstractLoadedObject end

struct LoadedStep <: AbstractLoadedObject
    state::AbstractState
    step::Int
    time::Float64

    function LoadedStep(state::AbstractState, step::Int, time::Float64)
        new(state, step, time)
    end
    function LoadedStep(input::Tuple)
        state, step, time = input
        new(state, step, time)
    end
end

struct LoadedSamplers <: AbstractLoadedObject

    n_samplers::Int
    info::Dict
    samplers::Vector

    function LoadedSamplers(n_samplers::Int, info::Dict, samplers::Vector)
        new(n_samplers, info, samplers)
    end


    function LoadedSamplers(input::Vector{<:AbstractSampler})

        samplers = input
        n_samplers = length(input)
        info = Dict{String, Vector}([])

        names = @. String(nameof(typeof(samplers)))
        unames = unique(names)
        for n in unames
            info[n] = []
        end
        

        for i in eachindex(samplers)
            sampler = samplers[i]
            name = String(nameof(typeof(sampler)))
            
            if name == "PointSampler"
                append!(info[name], sampler.sample_point_id)
            elseif name == "SectionSampler"
                coord = sampler.coord
                axis  = sampler.axis
                append!(info[name], [(axis=axis, coord=coord)])
            else
            end
        end

        new(n_samplers, info, samplers)

    end
end

function loadObj(input::Tuple)
    return LoadedStep(input)
end

function loadObj(input::AbstractCatalog)
    return input
end

function loadObj(input::AbstractCatalog, n_events::Int)
    catalog = input.catalog[1:n_events, :]

    return Catalog(catalog)
end

function loadObj(input::Vector{<:AbstractSampler})

    return LoadedSamplers(input)
end
function loadObj(input::AbstractSampler)

    return LoadedSamplers([input])
end

function loadData(input::String)
    data = load(input)["data"]
    return loadObj(data)
end

function loadData(input::String, n_events::Int)
    data = load(input)["data"]
    return loadObj(data, n_events)
end

function loadSSH(url, username, private_file, public_file)
    # url = ENV["elja_url"]
    # username = ENV["elja_user"]
    # private_file = ENV["elja_private"]
    # public_file = ENV["elja_pub"]
    sftp = SFTP(url, username, public_file, private_file)


    data  = load(SFTPClient.download(sftp, catalog_path))["data"]
    return loadObj(data)


end

