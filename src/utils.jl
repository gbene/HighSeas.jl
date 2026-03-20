"""
memcopy(A::AbstractGPUArray)
memcopy(A::AbstractArray, dev=0)

Copy an array to and from GPU/CPU depending on the output data type.

Overload this method with GPU specific Arrays using extensions

# Examples

A = rand(256, 256)
memcopy(A) # Copy A to the default device of an available backend
memcopy(A, 1) # Copy A to device 1 of the available backend
"""
memcopy(A::AbstractGPUArray{T,N}, dev_id::Int=0) where {T, N} = Array{T, N}(A)

function removenans(data::Vector)
      mask = any(!isnan, data; dims=2)
      return data[mask]
end


function removenans(data::Matrix)
      mask = any(!isnan, data; dims=2)
      return data[mask, :]
end

# function mycopy(x::T) where T

#       fnames = fieldnames(T)
#       v = Vector{Any}(undef, length(fnames))
#       for i in eachindex(fnames)
#             fname = fnames[i]
#             prop = deepcopy(getfield(x, fname))
#             prop_type = typeof(prop)
#             if prop_type <: Vector || prop_type <: Matrix
#                   prop = removenans(prop)
#             end
#             v[i] = prop
#       end
#       return T(v...)
# end


function ReadSheet(path_file::String; factor::Int=1)
      input_dict = Dict{String, Any}()
      open(path_file) do f
            lines = readlines(f)
            mask = @. !(isempty(lines) || contains(lines, "#"))
            lines = lines[mask]
            for line in lines
                  key, value = split(line,':')
                  if key in ["W","L","Wf","Lf","w","l","h","xi","yi","wi","li"]
                    value = parse(Float64, value)/factor
                  else
                    value = parse(Float64, value)
                  end
                  # display(key)
                  push!(input_dict,key=>value)
            end
      end
      input_dict
end

"""
CalcMinDt(fract, cs, cellsize)
CalcMinDt(input_dict::Dict)

Calculate minimum dt given the wave speed and cellsize
"""
function CalcMinDt(fract, cs, cellsize)

    frac     = 1/2^(fract-1)

    if frac <= 1.0
        mindt = frac/cs*2*cellsize;
    else
        mindt = 1/cs*2*cellsize;
    end
    return mindt
end
function CalcMinDt(input_dict::Dict)

      fract    = input_dict["fract"]
      cs       = input_dict["cs"]
      gridside = input_dict["cellsizex"]
      frac     = 1/2^(fract-1)

      if frac <= 1.0
            mindt = frac/cs*2*gridside;
      else
            mindt = 1/cs*2*gridside;
      end

      return mindt
end

"""
CheckLengthScales(material::AbstractMaterial, domain::AbstractDomain, σ::Float64)
CheckLengthScales(input_dict::Dict, material::AbstractMaterial, domain::AbstractDomain)


Check if grid and smallest dimension of the RW domain are within the imposed ratios:
    Lb/gridside > 3, Linf/gridside>10, L/Linf > 1
"""
function CheckLengthScales(material::AbstractMaterial, domain::AbstractDomain, σ::Float64)


      G = material.G
      ν = material.nu
      a = material.a
      b = material.b
      gridside = domain.grid.cell_sizex
      Dc = material.Dc
      L = domain.patch.l
      H = domain.patch.w

      Lb = G/(1-ν)*Dc/(b*σ)

      Linf = π/4*Lb*(b/(b-a))^2

      length_dim = min(H, L)
      ratio = length_dim/Linf

      if gridside*3 > Lb
            error("Cohesive zone may be poorly resolved, $(gridside*3) > $Lb")
      elseif gridside*10 > Linf
            error("Linf is poorly resolved,  $(gridside*10) > $Linf")
      elseif ratio < 1
            error("Linf is poorly resolved. H/Linf: $ratio")
      else
            return (Lb=Lb, Linf=Linf, ratio=ratio)
      end
end
function CheckLengthScales(input_dict::Dict, material::AbstractMaterial, domain::AbstractDomain)

    si = input_dict["si0"]

    res = CheckLengthScales(material, domain, si)
    return res
end


function RandomState(Nrows, Ncols)

    dx = rand(Ncols, Nrows)
    V  = rand(Ncols, Nrows)
    theta = rand(Ncols, Nrows)
    tau = rand(Ncols, Nrows)

    state = State(dx, V, theta, tau)
    return state
end


"""
Internal function used to set the backend used to perform the calculations.
Use the publicly available set_CPUbackend and set_GPUbackend to properly set the desired backend.
"""
function set_backend(backend::AbstractBackend)
    platform = backend.platform

    if platform in supported_platforms
        println(styled"{bold:$platform} will now be used")
        global_settings["backend"] = backend
    else
        error(styled"Platform {bold:$platform} is not supported")
    end
    return nothing
end

"""

get_backend()

Return the current backend used to perform the calculations.

# Example

```julia
get_backend()

>>> "CPUBackend"
```
"""
function get_backend()
    return global_settings["backend"]
end


"""
Get the list of supported platforms that can be used for running the simulations

# Example

```julia
get_available_platforms()

>>> ["CPU", "CUDA", ...]
```
"""
get_available_platforms() = display(supported_platforms)

"""
Get the list of supported GPU platforms that can be used for running the simulations

# Example

```julia
get_available_GPUplatforms()

>>> ["CUDA", ...]
```
"""
get_available_GPUplatforms() = display(supported_GPU_platforms)


"""
set_CPUbackend()

Set the backend to CPU

# Example

```julia
set_CPUbackend()

>>> CPU will now be used
```
"""
set_CPUbackend() = set_backend(CPUBackend())

"""
set_GPUbackend()

Set the backend to the available GPU. This is decided by the installed packages, i.e. if CUDA.jl is
installed (and imported) then CUDA will be used. Only import one JuliaGPU package per script.
**If no JuliaGPU package is loaded, this method will return a LoadError.**

See the available GPU platforms by running get_available_GPUplatforms().

# Example

```julia
set_GPUbackend()

>>> CUDA with CUDA.DeviceMemory will now be used
```

# Notes

Some GPU backends support UnifiedMemory (e.g. CUDA, ROCm, Metal). The default will always be DeviceMemory but
users can choose UnifiedMemory by running ```set_GPUbackend("unified")```
"""
function set_GPUbackend end


function make_outdir(start_time, output_dir)

    start_time = replace(start_time, ":"=>"_", "-"=>"_")
    device = get_backend().device

    outpath = "$(output_dir)/$(start_time)/$(device)"
    if ~isdir(outpath)
        mkpath(outpath)
    end

    return outpath

end
