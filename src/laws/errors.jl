

"""
DoubleError(experiment)

Create a DoubleError object used to calculate two errors between matrix A and B:

+ err (local error):  norm(A-B)/norm(A)
+ erra (absolute error): maximum(abs(A-B))/A[A==maximum(abs(A-B))]

err and erra are returned.

# Attribures

+ diff::AbstractArray{Float64}              # Matrix to store the difference A-B
+ absdiff::AbstractArray{Float64}           # Matrix to store the abs(A-B)
+ mask::AbstractArray{Bool}                 # Matrix to store where A==maximum(abs(A-B))

"""
struct DoubleError{M<:AbstractArray{Float64}, B<:AbstractArray{Int8}} <: AbstractErrorLaw

    diff::M
    absdiff::M
    mask::B

    function DoubleError(experiment::AbstractExperiment; gpu_id::Int=0)

        Nrows = experiment.domain.grid.n_elementsy
        Ncols = experiment.domain.grid.n_elementsx

        diff        = zeros(Nrows, Ncols)
        absdiff     = zeros(Nrows, Ncols)
        mask        = zeros(Int8, Nrows, Ncols)

        if typeof(get_backend()) <: AbstractGPUBackend
            diff        = memcopy(diff, gpu_id)
            absdiff     = memcopy(absdiff, gpu_id)
            mask        = memcopy(mask, gpu_id)
        end

        new{typeof(diff), typeof(mask)}(diff, absdiff, mask)
    end

    function DoubleError{M, B}(diff::M, absdiff::M, mask::B) where {M, B}
        new{M, B}(diff, absdiff, mask)
    end

end

function (calcDoubleError::DoubleError)(A, B)
    diff = calcDoubleError.diff
    absdiff = calcDoubleError.absdiff
    mask = calcDoubleError.mask

    @.. thread=true diff = A-B
    @.. thread=true absdiff = abs(diff)

    err = norm(diff)/norm(A)
    erra = maximum(absdiff)

    # This should make "indexing" GPU compatible code
    @.. thread=true mask = absdiff == erra # Using bool this allocates because thread=true does something to types. I think that the bool evaluation outputs Int so it gets converted to bool

    @.. thread=true absdiff = (erra / A) * mask # Using bool this allocates because it does type promotion to Int

    erra = sum(absdiff)

    return err, erra
end
