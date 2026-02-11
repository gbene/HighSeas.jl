module CudaExt

    using CUDA
    using HighSeas
    using StyledStrings

    function HighSeas.set_GPUbackend(mem::String="device")

        if mem == "device"
            backend = HighSeas.CUDABackend("GPU", "CUDA", CUDA.DeviceMemory)

        elseif mem == "unified"
            backend = HighSeas.CUDABackend("GPU", "CUDA", CUDA.UnifiedMemory)

        else
            error(styled"Memory {bold:$mem} type not recognized")
        end
        println(styled"CUDA with {bold:$(backend.memory)} will now be used")



        global_settings["backend"] = backend
        return nothing
    end


    function HighSeas.memcopy(A::AbstractArray{T, N}, dev_id::Int=0) where {T, N}
        mem = HighSeas.get_backend().memory

        prev_dev = device() # save current device

        device!(dev_id) # change to selected device

        A_cu = CUDA.CuArray{T, N, mem}(A) # move to memory

        device!(prev_dev) # change back to starting device

        return A_cu
    end

    function HighSeas.memcopy(A::AbstractArray{T, N}, mem::DataType, dev_id::Int=0) where {T, N}

        prev_dev = device() # save current device

        device!(dev_id) # change to selected device

        A_cu = CUDA.CuArray{T, N, mem}(A) # move to memory

        device!(prev_dev) # change back to starting device

        return A_cu
    end

end
