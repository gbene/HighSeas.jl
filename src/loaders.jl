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

function loadObj(input::Tuple)
    return LoadedStep(input)
end

function loadObj(input::AbstractCatalog)
    catalog = input.catalog
    mask = any(isnan, catalog; dims=2)
    catalog = catalog[.!vec(mask), :]


    return Catalog(catalog)
end

function loadObj(input::AbstractCatalog, n_events::Int)
    catalog = input.catalog[1:n_events, :]

    return Catalog(catalog)
end

function loadData(input::String)
    data = load(input)["data"]
    return loadObj(data)
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

