using SURD
using Documenter

DocMeta.setdocmeta!(SURD, :DocTestSetup, :(using SURD); recursive=true)

makedocs(;
    modules=[SURD],
    authors="Takeru Hashimoto",
    sitename="SURD.jl",
    format=Documenter.HTML(;
        canonical="https://tkrhsmt.github.io/SURD.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "API.md",
    ],
)

deploydocs(;
    repo="github.com/tkrhsmt/SURD.jl",
    devbranch="main",
)
