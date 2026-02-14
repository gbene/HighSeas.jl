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
    samplers_info::Vector
    samplers::Vector

    function LoadedSamplers(state::AbstractState, step::Int, time::Float64)
        new(state, step, time)
    end


    function LoadedSamplers(input::Vector{<:AbstractSampler})

        samplers = input
        n_samplers = length(input)
        samplers_info = Vector{NamedTuple}(undef, n_samplers)

        for i in eachindex(samplers_info)
            sampler = samplers[i]
            name = String(nameof(typeof(sampler)))
            if name == "PointSampler"
                point_id = sampler.sample_point_id
                samplers_info[i] = (name=name, point_id=point_id)
            elseif name == "SectionSampler"
                coord = sampler.coord
                axis  = sampler.axis
                samplers_info[i] = (name=name, axis=axis, coord=coord)
            else
            end
        end

        new(n_samplers, samplers_info, samplers)

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

