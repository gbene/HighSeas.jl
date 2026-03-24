using Documenter, HighSeas, CairoMakie

pages = ["Home" => "index.md",
         "Basic function and Objects" => "basic.md",
         "Materials" => "materials.md",
         "Grids, Geometries and Domains" => "geometries.md",
         "Laws" => "laws.md",
         "Algorithm objects" => "algorithms.md",
         "Detector objects" => "detectors.md",
         "Experiments" => "experiments.md",
         "Samplers" => "samplers.md",
         "Solving" => "solvers.md",
         "Saving and Reading" => "io.md",
         "Backends" => "backends.md",
         "Code design" => "structure.md"]

makedocs(sitename="HighSeas.jl", pages=pages, clean=true)


deploydocs(repo="github.com:gbene/HighSeas.jl.git")
