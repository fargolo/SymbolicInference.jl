using SymbolicInference
using Documenter

DocMeta.setdocmeta!(SymbolicInference, :DocTestSetup, :(using SymbolicInference); recursive=true)

makedocs(;
    modules=[SymbolicInference],
    authors="= <=> and contributors",
    repo=Documenter.Remotes.GitHub("fargolo", "SymbolicInference.jl"),
    sitename="SymbolicInference.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fargolo.github.io/SymbolicInference.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "RQA Moments Vignette" => "vignette_moments.md",
    ],
)

deploydocs(;
    repo="github.com/fargolo/SymbolicInference.jl",
    devbranch="main",
    push_preview=true,
)
