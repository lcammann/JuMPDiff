push!(LOAD_PATH, "../src/");
using Documenter, JuMPDiff

makedocs(sitename="JuMPDiff",
         authors = "Lucas Cammann",
         pages = ["Introduction" => "index.md"]
);