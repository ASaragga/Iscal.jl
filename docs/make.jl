using Iscalc
using Documenter

DocMeta.setdocmeta!(Iscalc, :DocTestSetup, :(using Iscalc); recursive=true)

makedocs(;
    modules=[Iscalc],
    authors="António Saragga Seabra",
    repo="https://github.com/ASaragga/Iscalc.jl/blob/{commit}{path}#{line}",
    sitename="Iscalc.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ASaragga.github.io/Iscalc.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ASaragga/Iscalc.jl",
)
