
"""
    StressFFT <: AbstractStressLaw

Analytical stress-slip relationship in the Fourier/wavenumber domain.

### Fields

+ `akdhat::AbstractArray{Float64}` -- akdhat values
+ `taur::AbstractArray{Float64}` -- Matrix used to store the calculated stress values from fft
+ `tau::AbstractArray{Float64}` -- Matrix used to store the final stress (tau = tau0+taur)
+ `F::AbstractArray{ComplexF64}` -- Matrix used to store the matrix mul result in the wavenumber domain
+ `dx_hat::AbstractArray{ComplexF64}` -- Matrix used to store the slip in the wavenumber domain
+ `p_rfft::AbstractFFTs.Plan` -- FFTW plan for forward FFT
+ `p_irfft::AbstractFFTs.Plan` -- FFTW plan for inverse FFT

### Notes

- When using GPUs, it is possible to decide where the masks reside using `gpu_id`
- We use rfft (https://juliamath.github.io/AbstractFFTs.jl/stable/api/#AbstractFFTs.rfft)
to be more memory efficient since we deal only with real numbers (i.e. slip can't be complex).
Also we don't specify the strategy for finding a plan in FFTW because cuFFT ignores it so it is useless.

### Examples

- `StressFFT(experiment::AbstractExperiment; gpu_id::Int=0)` -- Stress law for the current experiment
"""
struct StressFFT{M<:AbstractArray{Float64}, C<:AbstractArray{ComplexF64}} <: AbstractStressLaw

    akdhat::M
    taur::M
    F::C
    dx_hat::C
    p_rfft::AbstractFFTs.Plan
    p_irfft::AbstractFFTs.Plan

    tau::M

    function StressFFT(experiment::AbstractExperiment; gpu_id::Int=0)
        gridside = experiment.domain.grid.cell_sizex
        Nx = experiment.domain.grid.n_elementsx
        Ny = experiment.domain.grid.n_elementsy

        nu = experiment.material.nu
        G = experiment.material.G

        n_plan = div(Ny,2)+1
        Fs = 1/gridside;
        freqX = Fs*(-Nx/2:1:Nx/2-1)/Nx;
        freqY = Fs*(-Ny/2:1:Ny/2-1)/Ny;


        freqx = freqX'.*ones(Ny,Nx)
        freqy = freqY.*ones(Ny,Nx)



        kvx = ifftshift(2*pi*freqx)[1:n_plan,:];
        kvy = ifftshift(2*pi*freqy)[1:n_plan,:];

        antiplane = G;
        inplane = G/(1-nu);
        denom = sqrt.((kvx).^2 + (kvy).^2);

        akdhatx=   -inplane/2*(kvx).^2 ./denom;
        akdhaty=   -antiplane/2*(kvy).^2 ./denom;


        akdhatx[isnan.(akdhatx)] .= 0;
        akdhaty[isnan.(akdhaty)] .= 0;
        akdhat = (akdhatx + akdhaty);

        taur = zeros(Ny, Nx)
        F = zeros(ComplexF64, n_plan, Nx)

        dx_hat  = F.*akdhat

        p_rfft  = plan_rfft(taur)
        p_irfft = plan_irfft(dx_hat, Ny)

        tau = experiment.state.tau

        if typeof(get_backend()) <: AbstractGPUBackend
            akdhat      = memcopy(akdhat, gpu_id)
            taur        = memcopy(taur, gpu_id)
            F           = memcopy(F, gpu_id)
            dx_hat      = memcopy(dx_hat, gpu_id)
            p_rfft      = plan_rfft(taur)
            p_irfft     = plan_irfft(dx_hat, Ny)
        end

        new{typeof(tau), typeof(F)}(akdhat, taur, F, dx_hat, p_rfft, p_irfft, tau)

    end

    function StressFFT{M, C}(akdhat::M, taur::M, F::C, dx_hat::C, p_rfft, p_irfft, tau::M) where {M, C}
        new{M, C}(akdhat, taur, F, dx_hat, p_rfft, p_irfft, tau)
    end


end

function (stressFFT::StressFFT)(dx, tau0)
    F = stressFFT.F
    p_fft = stressFFT.p_rfft
    p_ifft = stressFFT.p_irfft
    tau = stressFFT.tau

    dx_hat = stressFFT.dx_hat
    akdhat = stressFFT.akdhat
    taur = stressFFT.taur
    mul!(F, p_fft, dx) #Matrix multiplication using BLAS (or CUBLAS). This is the fft, i.e. F = rfft(dx)
    @.. thread=true dx_hat = akdhat*F
    mul!(taur, p_ifft, dx_hat) # This is the ifft i.e. F = irfft(dx_hat, Ny)

    @.. thread=true tau = tau0+taur;
    return tau

end
