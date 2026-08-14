
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

- `catalog::Matrix{Float64}` -- Matrix representation of the catalog
- `max_events::Int` -- Maximum number of catalog events (i.e. number of rows)
- `t::SubArray` -- Time of rupter (in seconds)
- `interevent_time::SubArray` -- Interevent time between (in seconds)
- `Moment::SubArray` -- Moment released (in N⋅m)
- `mag::SubArray` -- Magnitude using Hanks and Kanamori 1979
- `Area::SubArray` -- Rupture area (m²)
- `MeanSlip::SubArray` -- Mean slip in the rupture area (m)
- `MeanStress::SubArray` -- Mean stress drop in the rupture area (Pa)
- `hypo_x::SubArray` -- x coordinate of the event (m)
- `hypo_y::SubArray` -- y coordinate of the event (m)
- `n_events::Int` -- Numbe of recorded events



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

struct EventLog <: AbstractEvent
    event_n::Int
    event_time::Float64
    event_mag::Float64
    starttimestamp::String
    endtimestamp::String

end

struct SimulationLog <: AbstractLog
    exp_starttime:: String
    input_dict::Dict
    length_scales::String
    events::Vector{<:EventLog}
    nevents::Int


    function SimulationLog(exp_starttime::String, input_dict::Dict, length_scales::String, events::Vector{<:EventLog}, n_events::Int)
        new(exp_starttime, input_dict, length_scales, events, n_events)
    end


    function SimulationLog(log_file_path::String)
        lines = readlines(log_file_path)

        exp_starttime = split(lines[1], ": ")[2]
        input_values = split.(lines[3:29],": ")
        input_dict = Dict{String, Float64}(k=>parse(Float64, v) for (k,v) in input_values)

        length_scales = lines[31]

        n_lines_events = length(lines)
        if isodd(n_lines_events) # if it's odd it means that the simulation stopped while the event was going on
            n_lines_events -= 1
        end

        event_lines = 33:2:n_lines_events
        events = Vector{EventLog}(undef, length(event_lines))

        for (a,i) in enumerate(event_lines)
            vals = split(lines[i], " ")
            vals2 = split(lines[i+1], " ")
            starttimestamp = vals[1]
            event_n = parse(Int, vals[3])
            event_time=parse(Float64, vals[end])
            event_mag=parse(Float64, vals2[end])
            endtimestamp = vals[1]
            event = EventLog(event_n, event_time, event_mag, starttimestamp, endtimestamp)

            events[a] = event


        end
        n_events = length(events)
        new(exp_starttime, input_dict, length_scales, events, n_events)
    end



end


"""
MailNotifier

Setup a mail notification system. 

An email will be sent to the desired recipient(s) at the start and end of a simulation


### Fields

- `mailto::Vector{String}` -- Recipent(s) that will receive the email
- `mailfrom::String` -- Sender of the email
- `url::String` -- SMTP url
- `send_opts::SendOptions` -- Send options struct


### Notes

For security reasons the sender email and password must be set as environment variables. When constructing the object, pass the names of the env variables for the email username and password. 


### Examples

- `MailNotifier(["mailto@example.com"], "name_of_env_variable_user", "name_of_env_variable_password", "smtps://smtp.gmail.com:465")` -- Mail notifier that uses gmail, sending the email to a receiver 
- `MailNotifier("name_of_env_variable_user", "name_of_env_variable_password", "smtps://smtp.gmail.com:465")` -- Mail notifier that uses gmail, sending the email to yourself
- `MailNotifier("name_of_env_variable_user", "name_of_env_variable_password", "smtps://smtp.gmail.com:465"; test=true)` -- Test the mail notifier by sending a test email (default is false)
- `MailNotifier("name_of_env_variable_user", "name_of_env_variable_password", "smtps://smtp.gmail.com:465"; isSSL=false)` -- Don't use SSL certificates (default is true)

"""
struct MailNotifier

    mailto::Vector{String}
    mailfrom::String
    url::String
    
    send_opts::SendOptions

    function MailNotifier(mailto::Vector{String}, env_mail_user::String, env_mail_pass::String, url::String; test::Bool=false, isSSL=true)

        mail_user=ENV[env_mail_user]
        mail_pass=ENV[env_mail_pass]
        opt = SendOptions(isSSL = isSSL, username = mail_user, passwd = mail_pass)

        notifier = new(mailto, mail_user, url, opt)

        if test 
            send(notifier, "Test email", "This is a test email to check if the notifier works")
        end

        return notifier
    end

    function MailNotifier(env_mail_user::String, env_mail_pass::String, url::String; test::Bool=false, isSSL=true)

        mail_user=ENV[env_mail_user]
        mail_pass=ENV[env_mail_pass]
        opt = SendOptions(isSSL = isSSL, username = mail_user, passwd = mail_pass)

        mailto = [mail_user]

        notifier = new(mailto, mail_user, url, opt)

        if test 
            send(notifier, "Test email", "This is a test email to check if the notifier works")
        end
        
        return notifier
    end

end


"""
send(notifier, "This is the subject", "This is a message")
send(notifier, "This is the subject", "This is a message", ["path/to/attachment1.png", "path/to/attachment2.txt"])

Send an email using a MailNotifier object.

This function sends an email using the options of a MailNotifier object. Attachments can be added when providing a vector of paths.

### Arguments

- `notifier::MailNotifier` -- MailNotifier object
- `subject::String` -- Subject of the email
- `message::String` -- Message of the email
- `attachments::Vector{String}` -- Attachments of the email
"""
function send(notifier::MailNotifier, subject::String, message::String)

    
    body = get_body(notifier.mailto, notifier.mailfrom, subject, get_mime_msg(message))

    resp = SMTPClient.send(notifier.url, notifier.mailto, notifier.mailfrom, body, notifier.send_opts)

    return resp

end

function send(notifier::MailNotifier, subject::String, message::String, attachments::Vector{String})

    
    body = get_body(notifier.mailto, notifier.mailfrom, subject, get_mime_msg(message); attachments)

    resp = SMTPClient.send(notifier.url, notifier.mailto, notifier.mailfrom, body, notifier.send_opts)

    return resp

end
