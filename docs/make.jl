using MonolithicFEMVLFS
using Documenter

DocMeta.setdocmeta!(MonolithicFEMVLFS, :DocTestSetup, :(using MonolithicFEMVLFS); recursive=true)

makedocs(;
    modules=[MonolithicFEMVLFS],
    authors="Oriol Colomes <oriol.colomes@gmail.com>",
    repo="https://github.com/oriolcg/MonolithicFEMVLFS.jl/blob/{commit}{path}#{line}",
    sitename="MonolithicFEMVLFS.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://oriolcg.github.io/MonolithicFEMVLFS.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/oriolcg/MonolithicFEMVLFS.jl",
    devbranch="main",
)
