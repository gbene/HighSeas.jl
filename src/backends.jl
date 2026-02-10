abstract type AbstractBackend end
abstract type AbstractGPUBackend<:AbstractBackend end


struct CPUBackend <: AbstractBackend

    platform::String
    device::String

    function CPUBackend()
        new("CPU", "CPU")
    end

end

struct CUDABackend <: AbstractGPUBackend

    platform::String
    device::String
    memory::DataType

end

struct MetalBackend <: AbstractGPUBackend

    platform::String
    device::String
    memory::DataType

end
