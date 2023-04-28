using AngularSpectrumMethod
using Documenter

DocMeta.setdocmeta!(AngularSpectrumMethod, :DocTestSetup, :(using AngularSpectrumMethod); recursive=true)

makedocs(;
    modules=[AngularSpectrumMethod],
    authors="Shuhei Yoshida <yshuhei@ele.kindai.ac.jp> and contributors",
    repo="https://github.com/syoshida1983/AngularSpectrumMethod.jl/blob/{commit}{path}#{line}",
    sitename="AngularSpectrumMethod.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://syoshida1983.github.io/AngularSpectrumMethod.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/syoshida1983/AngularSpectrumMethod.jl",
    devbranch="master",
)
