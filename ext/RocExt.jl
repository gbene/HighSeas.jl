module RocExt

    using AMDGPU
    using HighSeas
    using StyledStrings

    function HighSeas.set_GPUbackend(mem::String="device")

        if mem == "device"
            backend = HighSeas.ROCmBackend("GPU", "ROCm", mem, AMDGPU.Runtime.Mem.HIPBuffer)

        elseif mem == "unified"
            # backend = HighSeas.MetalBackend("GPU", "METAL", mem, Metal.SharedStorage)
            error(styled"Memory {bold:$mem} not supported")


        else
            error(styled"Memory {bold:$mem} type not recognized")
        end
        println(styled"AMDGPU with {bold:$(backend.memory)} will now be used")



        global_settings["backend"] = backend
        return nothing
    end

    # AMDGPU has a weird API. device and device_id are distinct functions and indexing is 1 based.
    # To be uniform with the others API, the input dev_id must be 0 based and so we adapt the code
    # correspondingly!
    function HighSeas.memcopy(A::AbstractArray{T, N}, dev_id::Int=0) where {T, N}
        mem = HighSeas.get_backend().memory

        dev_id +=1

        prev_dev = device() # save current device
        if device_id() != dev_id
            println(styled"Switching to device {bold:$(dev_id-1)}")

            device_id!(dev_id) # change to selected device

            A_roc = ROCArray{T, N, mem}(A) # move to memory

            println(styled"Switching back to orginal device {bold:$(device_id(prev_dev))}")
            device!(prev_dev) # change back to starting device
            return A_roc
        else
            A_roc = ROCArray{T, N, mem}(A) # move to memory
            return A_roc
        end
    end

end
