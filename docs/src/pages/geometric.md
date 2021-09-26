```@meta
DocTestSetup = quote
    using JuMPDiff, Ipopt
end
```
```@setup plot
using Plots; pyplot()
Plots.reset_defaults()
```

# Geometric example
Let us consider the function $f(x,y) = (2x^2 +3y^2) +4xy$ which we wish to minimize, with all candidate solutions required to lie on the unit circle $g(x,y) = x^2 + y^2 - 1 = 0$. We plot the objective function and the constraint for a better understanding of what we are working with. 

```@example plot
using Plots
x = -2:0.1:2;
y = -2:0.1:2;
f(x,y) = 2*x.^2 + 3*y.^2 + 4*x.*y    #Our objective function
surface(x,y,f);                      #Surface plot of objective function
x = -1:0.1:1
y = sqrt.(1 .- x.^2)
plot!(x,y,f(x,y),c=:cyan,legend = false)
plot!(x,-y,f(x,-y),c=:cyan, legend = false)

```
---


Our next step is to analytically find and analyze the solutions of this problem with the aid of JuMPDiff functions. Afterwards, we will compare these results to ones obtained numerically with the Ipopt solver. 
## Analytical solution 
For the analytical solution of this problem we will use the method of Lagrange multipliers. 
```jldoctest examplegeo
julia> using JuMPDiff, Ipopt

julia> model = Model(Ipopt.Optimizer); 

julia> @variable(model,x, start = 1);

julia> @variable(model,y, start = 1);

julia> @variable(model, λ); #Lagrange multiplier

julia> @constraint(model,circle, x^2 + y^2 - 1 == 0);

julia> @expression(model,g,x^2 + y^2 -1); #Left hand side of the constraint

julia> @objective(model, Min, (2x^2 +3y^2) +4x*y); 
```
The Lagrangian of this problem is defined as 
```math 
\begin{aligned}
\mathcal{L}(x,y,\lambda) = f(x,y) - \lambda g(x,y)
\end{aligned}
```
and we seek to find $x$, $y$ and $\lambda$ s.t. 
```math 
\begin{aligned}
\nabla_{x,y,\lambda} \mathcal{L}(x,y,\lambda) = \nabla_{x,y}f(x,y) - \lambda \nabla_{x,y} g(x,y) = 0
\end{aligned} 
```
The partial derivatives of $\mathcal{L}$ with respect to $x$ and $y$ can be easily constructed with the `con_der` and `obj_der` functions, the partial derivative w.r.t. $\lambda$ is merely $g$ itself. 
```jldoctest examplegeo
julia> ∇L = obj_der(model,[x,y]) - λ*con_der(model,circle,[x,y])
2-element Array{GenericQuadExpr{Float64,VariableRef},1}:
 -2 λ*x + 4 x + 4 y
 -2 λ*y + 6 y + 4 x

julia> push!(∇L, g)
3-element Array{GenericQuadExpr{Float64,VariableRef},1}:
 -2 λ*x + 4 x + 4 y
 -2 λ*y + 6 y + 4 x
 x² + y² - 1
```
 
!!! note 
    It is also possible to first construct the Lagrangian and then only use `obj_der` to compute its gradient. In the above example we were however required to assemble the gradient from the derivatives of the constituting expressions, since the constraint is quadratic and can not be multiplied with λ without making it a `NLexpression`.

For the sake of brevity the analytical solution procedure for solving $∇\mathcal{L} = 0$ shall rather be briefly explained than shown in great detail. One of the ways of solving the set of equations is to first eliminate $λ$ from the first two equations. By then isolating $x$ or $y$ in the third equation and inserting the resulting expression in equations one and two we can find numerical values for either one of these variables. Doing so we find
```math 
\begin{aligned}
y_1 = x_1 = \pm \sqrt{\frac{1}{2} - \frac{1}{2 \sqrt{17}}} \approx \pm 0.61541 \\
y_2 = x_2 = \pm \sqrt{\frac{1}{2} + \frac{1}{2 \sqrt{17}}} \approx \pm  0.78821 \\
\end{aligned}
```
where all combinations of $x$ and $y$ with different indices are valid critical points. To learn more about the nature of these critical points we can evaluate the bordered Hessian $\bar{H}$, which for the current case is defined as 
```math 
\begin{aligned}
\bar{H}(x,y) =\begin{bmatrix}0&{\frac  {\partial g}{\partial x}}&{\frac  {\partial g}{\partial y}}\\[1.5ex]{\frac  {\partial g}{\partial x}}&{\frac  {\partial ^{2}\mathcal{L}}{\partial x^{2}}}&{\frac  {\partial ^{2}\mathcal{L}}{\partial x\partial y}}\\[1.5ex]{\frac  {\partial g}{\partial y}}&{\frac  {\partial ^{2}\mathcal{L}}{\partial y\partial x}}&{\frac  {\partial ^{2}\mathcal{L}}{\partial y^{2}}}\\\end{bmatrix}
\end{aligned}
```
We can easily construct the bordered Hessian by using the `model_hess` function and the `con_der` function. In doing so we can also make use of the fact that JuMPDiff offers the possibility to have derivatives evaluated at specified points. 
```jldoctest examplegeo
julia> H = Array{GenericAffExpr{Float64,VariableRef}}(undef,(3,3)); #Initialize array

julia> H[1,1] = 0;

julia> H[1,2:3] = con_der(model,circle,[x,y],[-0.7882, 0.6125]); #Evaluated at given point

julia> H                            
3×3 Array{GenericAffExpr{Float64,VariableRef},2}:
    0       -1.5764     1.225
 #undef  #undef      #undef
 #undef  #undef      #undef

julia> H[2:3,1] = con_der(model,circle,[x,y],[-0.7882, 0.6125]);

julia> H 
3×3 Array{GenericAffExpr{Float64,VariableRef},2}:
 0           -1.5764     1.225
 -1.5764  #undef      #undef
 1.225    #undef      #undef

julia> H[2:3,2:3] = model_hess(model,[x,y]);

julia> H
3×3 Array{GenericAffExpr{Float64,VariableRef},2}:
 0        -1.5764  1.225
 -1.5764  4        4
 1.225    4        6
```
The determinant of $\bar{H}(-0.7882,0.6125)$ turns out to be negative, and we can therefore conclude that this point is a local minimizer. Taking into account the curvature of the feasible set in the above figure we can even be certain that this is actually the global minimizer of our problem. 


!!! note 
    In fact, the point (0.7882,-0.6125) equally minimizes our problem due to its symmetry.

## Numerical Solution 
Clearly, the above example was artificial since we can find the minimum way easier numerically.
```jldoctest examplegeo
julia> set_silent(model) 
true 

julia> optimize!(model);

julia> value(x)
0.788205437944948

julia> value(y)
-0.6154122094937813
```
The solver finds one of the minima which we have previously identified analytically. To find the second minima, we would have to supply a different initial guess to `Ipopt`.

