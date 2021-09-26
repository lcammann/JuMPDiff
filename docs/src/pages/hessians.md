```@meta
DocTestSetup = quote
    using JuMPDiff
end
```

# Hessians 

## Basics 
The Hessian of a scalar-valued differentiable function $f$ is a square matrix of its second order partial derivatives. Hence, at a point $p = (x_1,...,x_n)$ for $f\colon \mathbb {R} ^{n}\to \mathbb {R}$ the Hessian $\mathbf{H}_f\colon \mathbb {R} ^{n}\to \mathbb {R}^{n \times n}$ is defined as : 
```math
\begin{aligned} 
\mathbf {H} _{f}(p) = {\begin{bmatrix}{\dfrac {\partial ^{2}f}{\partial x_{1}^{2}}(p)}&{\dfrac {\partial ^{2}f}{\partial x_{1}\,\partial x_{2}}(p)}&\cdots &{\dfrac {\partial ^{2}f}{\partial x_{1}\,\partial x_{n}}(p)}\\[2.2ex]{\dfrac {\partial ^{2}f}{\partial x_{2}\,\partial x_{1}}(p)}&{\dfrac {\partial ^{2}f}{\partial x_{2}^{2}}(p)}&\cdots &{\dfrac {\partial ^{2}f}{\partial x_{2}\,\partial x_{n}}(p)}\\[2.2ex]\vdots &\vdots &\ddots &\vdots \\[2.2ex]{\dfrac {\partial ^{2}f}{\partial x_{n}\,\partial x_{1}}(p)}&{\dfrac {\partial ^{2}f}{\partial x_{n}\,\partial x_{2}}(p)}&\cdots &{\dfrac {\partial ^{2}f}{\partial x_{n}^{2}}(p)}\end{bmatrix}}
\end{aligned}
```

JuMPDiff provides a function with different methods to compute the Hessian of the objective function of a JuMP model. 

## Overview 
Different dispatches of the function `model_hess` are available to compute the Hessian of the objective function registered in the JuMP model `model`. The basic functionalities are: 
```@docs 
model_hess
```
Results retrieved by these methods can also be registered in the model by calling the function `model_hess!`. 
```@docs
model_hess!
```
## Examples
The following shall illustrate the use of these methods. Note that the calls of `model_hess` below are independant of each and shall merely illustrate the different available methods. For the sake of the example, let us consider an arbitray objective function with four variables: 
```jldoctest examplehess
julia> model = Model();

julia> @variable(model,w);

julia> @variable(model,x);

julia> @variable(model,y);

julia> @variable(model,z);

julia> @objective(model, Max, 1*w + 2*x^2 + 3*x*y + 4*y*z)
2 x² + 3 y*x + 4 z*y + w

julia> model_hess(model)
4×4 Array{Float64,2}:
 0.0  0.0  0.0  0.0
 0.0  4.0  3.0  0.0
 0.0  3.0  0.0  4.0
 0.0  0.0  4.0  0.0

julia> model_hess(model, [x,y])
2×2 Array{Float64,2}:
 4.0  3.0
 3.0  0.0

julia> model_hess(model, [y,x])
2×2 Array{Float64,2}:
 0.0  3.0
 3.0  4.0
```
A few things appear noteworthy in the above examples. The first is the fact that JuMP prints the terms of the objective function in a different order than specified by the user. This is due to the fact that JuMP internally represents equations as combinations of pre-defined functions and sets, defined in the abstract data structure `MathOptInterface` (see also [^1]). Often times this entails equations to be rearranged.\
Secondly, we can see that the default variable order of the Hessian is the order in which the variables are registered in the model. In the above example this means that $x_1 = w$, $x_2 = x$ etc. Similarly, if the user supplies an array of variables as second function argument, the Hessian is evaluated according to the order in that array. For example, in the second call of `model_hess` above the variable array $\mathbf{x} = [x, y]$ is supplied, in which case the Hessian is computed as 
```math
\begin{aligned} 
\mathbf {H} _{f} = {\begin{bmatrix}{\dfrac {\partial ^{2}f}{\partial x^{2}}}&{\dfrac {\partial ^{2}f}{\partial x\,\partial y}} \\[2.2ex] 
{\dfrac {\partial ^{2}f}{\partial y\,\partial x}} & {\dfrac {\partial ^{2}f}{\partial y^{2}}}
\end{bmatrix}}
\end{aligned}
```
while for the variable array $\mathbf{x} = [y, x]$ supplied in the last function call the Hessian is [^2]
```math
\begin{aligned} 
\mathbf {H} _{f} = {\begin{bmatrix}{\dfrac {\partial ^{2}f}{\partial y^{2}}}&{\dfrac {\partial ^{2}f}{\partial y\,\partial x}} \\[2.2ex]
{\dfrac {\partial ^{2}f}{\partial x\,\partial y}} & {\dfrac {\partial ^{2}f}{\partial x^{2}}}
\end{bmatrix}}
\end{aligned}
```
Another fact to note is that the resulting array is of the type `Float64`, whereas for other functions included in JuMPDiff, such as `model_jac` or `con_der`, it is `GenericAffExpr{Float64,VariableRef}` by default. This can be explained by the fact that `JuMPDiff` only handles terms up to quadratic order, for which the second derivate is a scalar in any case. Other functions included in this package however compute first derivatives and therefore require a more flexible data type. 

The `model_hess!` function extends the previously described functionalities by registering the retrieved results within the JuMP model itself. 
```jldoctest examplehess
julia> model_hess!(model)
4×4 Array{Float64,2}:
 0.0  0.0  0.0  0.0
 0.0  4.0  3.0  0.0
 0.0  3.0  0.0  4.0
 0.0  0.0  4.0  0.0

julia> model[:Hessian]
4×4 Array{Float64,2}:
 0.0  0.0  0.0  0.0
 0.0  4.0  3.0  0.0
 0.0  3.0  0.0  4.0
 0.0  0.0  4.0  0.0

julia> model_hess!(model,[z,w],:Hessian_z_w)
2×2 Array{Float64,2}:
 0.0  0.0
 0.0  0.0

julia> model[:Hessian_z_w]
2×2 Array{Float64,2}:
 0.0  0.0
 0.0  0.0
```
The default name for the registered Hessian is `Hessian`. The result can then be called by querying `model[:Hessian]` (where model is the name of the JuMP model). By specifying a different name in the function arguments we can register the result under a different name, and make it callable under that name. This is what we did in the above example, since the second call of `model_hess!` would have otherwise attempted to re-register the default name `Hessian`, resulting in the following error: 
```jldoctest examplehess
julia> model_hess!(model,[z,w])
ERROR: You are trying to register the name Hessian in the model, which is already occupied. If you have not specified the name Hessian, then the function you have called assumes Hessian to be the default name. To overcome this issue specify a different name in the function arguments or unregister the name Hessian from the model.
[...]
```
A slightly different error would have been thrown if if were attempted to register the user-defined name twice:
```jldoctest examplehess
julia> model_hess!(model,[z,w], :Hessian_z_w)
ERROR: You are trying to register a name which is already occupied in the model. Choose a different one or unregister the existing.
[...]
```



[^1]: [Legat, Benoit and Dowson, Oscar and Garcia, Joaquim and Lubin, Miles     MathOptInterface: a data structure for mathematical optimization problems. _INFORMS Journal on Computing_, 2020 ](https://arxiv.org/abs/2002.03447)
[^2]: Note that this only influences the entries on the main diagonal, since the mixed partial derivates are equal by Schwarz's theorem