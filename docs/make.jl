using ScalarWaveFVM
using Documenter

DocMeta.setdocmeta!(ScalarWaveFVM, :DocTestSetup, :(using ScalarWaveFVM); recursive=true)

makedocs(;
    modules=[ScalarWaveFVM],
    repo=Documenter.Remotes.GitHub("svretina", "ScalarWaveFVM.jl"),
    authors="Stamatis Vretinaris",
    sitename="ScalarWaveFVM.jl",
    format=Documenter.HTML(;
        canonical="https://svretina.github.io/ScalarWaveFVM.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/svretina/ScalarWaveFVM.jl",
    devbranch="master",
)
