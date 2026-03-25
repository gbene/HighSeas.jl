![icon](docs/src/assets/banner.svg)


Ahoy, matey! Welcome to HighSeas, a Julia package focused on accelerating Sequences of Earthquakes and Aseismic Slips simulations using the Spectral Boundary Integral Method on GPUs. 

# How to install
Open the REPL and type either 

```using Pkg; Pkg.add(url="git@github.com:gbene/HighSeas.jl.git")``` 

or 

```] add git@github.com:gbene/HighSeas.jl.git``` 

Users are encouraged to use environments but it is not necessary. 


# How to use

Look at the examples folder for a complete example on how to use the library. Crude docs will eventually come. 

To run the examples, cd in the examples/CPU or GPU_CUDA directory and do either

```bash
julia --project=. --threads auto BP4QD.jl
```
to run on CPU

or
```bash
julia --project=. --threads auto BP4QD_GPU.jl
```
to run on GPU (CUDA)

Before trying both it is suggested to look at the code to understand what is being saved.

Then to plot the results, cd in examples/Plots and run an active REPL.

```bash
julia --project=. --threads auto 
```
then

```julia
include("BP4QD_plot.jl") 
```

To display the figure run

```julia
eventfig # or the other var names 
```

You can also find a MATLAB script version in `examples/resources/MATLAB`

# References

If you use this code please cite our work!

Benedetti, G., Heimisson, E. R. (2026, in preparation). _HighSeas.jl Accelerating rate and state simulations with GPUs_

# Acknowledgments

## Funding
This work was supported by Rannís IRF grant no. 2410397-051

## Libraries
We use many libraries of the Julia ecosystem and we thank all of them! These are the main ones that made this project possible

+ JuliaGPU: [https://juliagpu.org/](https://juliagpu.org/)
+ FFTW: 
    + [https://www.fftw.org/](https://www.fftw.org/) 
    + [https://github.com/JuliaMath/AbstractFFTs.jl](https://github.com/JuliaMath/AbstractFFTs.jl)
    + [https://github.com/JuliaMath/FFTW.jl](https://github.com/JuliaMath/FFTW.jl])
+ FastBroadcast.jl: [https://github.com/YingboMa/FastBroadcast.jl](https://github.com/YingboMa/FastBroadcast.jl)

## People
+ Elías Rafn Heimisson for desining the solver and providing constant support:  [https://github.com/eliasrh](https://github.com/eliasrh)
+ Ylse Anna de Vries for the name of the package: [https://github.com/ylseanna](https://github.com/ylseanna)