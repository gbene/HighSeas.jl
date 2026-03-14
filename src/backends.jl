


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
    memtype::String
    memory::DataType

end

struct ROCmBackend <: AbstractGPUBackend

    platform::String
    device::String
    memtype::String
    memory::DataType

end

struct MetalBackend <: AbstractGPUBackend

    platform::String
    device::String
    memtype::String
    memory::DataType

end
