push!(LOAD_PATH, "../src/");
using Documenter, JuMPDiff

makedocs(sitename="JuMPDiff",
         format = Documenter.HTML(
                            prettyurls = get(ENV, "CI", nothing) == "true"),
         authors = "Lucas Cammann",
         doctest = true,
         pages = ["Introduction" => "index.md",
                  "Gradients" => "pages/gradients.md",
                  "Jacobians" => "pages/jacobians.md",
                  "Hessians" => "pages/hessians.md",
                  "Geometric example" => "pages/geometric.md"]
);