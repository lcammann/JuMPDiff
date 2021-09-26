# JuMP-Diff

This package is intended to provide user friendly functions to calculate symbolic derivatives of JuMP models. It furthermore provides different case study files, in which the functions of the package are used to solve bilevel optimization problems. 

## Documentation 
The documentation provides an overview of the available functions, as well as examples on how they can be used. To view the documentation, you must have the Documenter package and build the documentation. This can be done as follows: 

1- Make sure you have the Documenter package. Press `]` in the julia REPL to start the Pkg REPL
```julia
(@v1.5) pkg> add Documenter
```

2- Run the `make.jl` file in the `docs` folder
```julia
julia> include("/path/to/JuMPDiff/docs/make.jl")
```

3- The html pages will be available on `/path/to/JuMPDiff/docs/build/`. Opening `index.html` will allow to navigate through the documentation.

## Case studies
Different case studies are available showing how `JuMPDiff` functions can be useful in solving different optimization problems, including Bilevel optimization problems. The first is inspired by the paper by [Mahadevan et al. (2002)](https://www.cell.com/biophysj/fulltext/S0006-3495(02)73903-9), which introduces the concept of dynamic flux balance analysis (DFBA) to model the batch growth of E. Coli on glucose. In short, the proposed approach allows for solving the dynamic reactor equations under consideration of the metabolic network within the microorganisms, which itself can mathematically be represented as an optimization problem. A comparison of different solution methods to this problem, including those using `JuMPDiff`, can be found in `Case study/DFBA-Functions.jl`.

![Batch problem](../master/docs/images/Batch.JPG)

Building on the initial problem, a second case study considers the same microbial growth in a CSTR. The problem statement is altered to be an optimization problem in which it is sought to maximize the overall throughput of biomass. Considering the overall internal metabolic network model makes this a bilevel optimization problem. The implementation of this problem can be found in `Case study/Bilevel-Functions.jl`. 

![CSTR problem](../master/docs/images/CSTR.JPG)

In a last case study it is sought to control the previously found optimal, time-dependant profile of the biomass concentration via the dilution rate by using a Model Predictive Controller (MPC). The corresponding functions can be found in `Case study/MPC-Functions.jl`. It is recommended to run and compare the respective files through the `Main.jl` file, as this loads the relevant dependencies and allows to conveniently alter the used process parameters.

