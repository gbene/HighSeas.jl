using Documenter, HighSeas, CairoMakie

pages = ["Home" => "index.md"]

makedocs(sitename="HighSeas.jl", pages=pages, clean=true)


deploydocs(repo="github.com:gbene/HighSeas.jl.git",devbranch="dev")
