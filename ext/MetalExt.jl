module MetalExt

    using Metal
    using HighSeas
    using StyledStrings

    # Metal PrivateStorage = CUDA DeviceMemory, Metal SharedStorage = CUDA UnifiedMemory
    function HighSeas.set_GPUbackend(mem::String="device")

        if mem == "device"
            backend = HighSeas.MetalBackend("GPU", "METAL", mem, Metal.PrivateStorage)

        elseif mem == "unified"
            backend = HighSeas.MetalBackend("GPU", "METAL", mem, Metal.SharedStorage)

        else
            error(styled"Memory {bold:$mem} type not recognized")
        end
        println(styled"METAL with {bold:$(backend.memory)} will now be used")



        global_settings["backend"] = backend
        return nothing
    end


    function HighSeas.memcopy(A::AbstractArray{T, N}, dev_id::Int=0) where {T, N}
        mem = HighSeas.get_backend().memory

        prev_dev = device() # save current device
        device!(dev_id) # change to selected device


        A_mtl = Metal.MtlArray{T, N, mem}(A)

        device!(prev_dev) # change to selected device

        return A_mtl
    end

end
