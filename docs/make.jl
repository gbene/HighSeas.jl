using Documenter, HighSeas, CairoMakie

pages = ["Home" => "index.md",
         "Grids, Geometries and Domains" => "geometries.md",
         "Algorithm objects" => "algorithms.md",
         "Detector objects" => "detectors.md",
         "Experiments" => "experiments.md",
         "Code design" => "structure.md"]

makedocs(sitename="HighSeas.jl", pages=pages, clean=true)


deploydocs(repo="github.com:gbene/HighSeas.jl.git",devbranch="dev")
