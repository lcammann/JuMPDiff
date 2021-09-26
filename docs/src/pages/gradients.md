```@meta
DocTestSetup = quote
    using JuMPDiff
end
```

# Gradients 

## Basics
The gradient of a scalar-valued differentiable function $f$ is a vector $\nabla f$ whose components are the partial derivatives of $f$. Hence, at a point $p = (x_1, ..., x_n)$ for $f\colon \mathbb {R} ^{n}\to \mathbb {R}$ the gradient $\nabla f\colon \mathbb {R} ^{n}\to \mathbb {R}^{n}$ is defined as : 
```math
\begin{aligned}
\nabla f(p)= {\begin{bmatrix} {\frac {\partial f}{\partial x_{1}}}(p)\\ \vdots \\{\frac {\partial f}{\partial x_{n}}}(p)\end{bmatrix}}
\end{aligned}
```
JuMPDiff provides different functions to calculate the derivative of the objective function and constraints of a JuMP model. 
## Objective function
### Overview
Different dispatches of the function `obj_der` are available to compute the gradient of the objective function registered in the JuMP-model `model`. The basic functionalities are the following:
```@docs 
obj_der
```
Results retrieved by these methods can also be registered in the model by calling the function `obj_der!`. 
```@docs
obj_der!
```

### Examples
The following shall illustrate the use of these methods: 

```jldoctest
julia> model = Model();

julia> @variable(model,x);

julia> @variable(model,y);

julia> @objective(model, Max, x*y - x^2 - y^2);

julia> obj_der(model)
2-element Array{GenericAffExpr{Float64,VariableRef},1}:
 y - 2 x
 x - 2 y

julia> obj_der(model,[x])
1-element Array{GenericAffExpr{Float64,VariableRef},1}:
 y - 2 x

julia> obj_der(model,[2,1])
2-element Array{GenericAffExpr{Float64,VariableRef},1}:
 -3
 0

julia> obj_der(model,[y],[3])
1-element Array{GenericAffExpr{Float64,VariableRef},1}:
 x - 6
```
Note that the array of variable values is evaluated according to the order in which the variables are registered in the model, i.e. in the third funtion call of `obj_der` the value 2 is associated with the variable `x` and the value 1 with the variable `y`. This order can be queried by the user through the command `index(var)`, in which `var` is the variable of which the load order is unknown. The returned index number corresponds to the belonging index in the `vals` array. Alternatively, the user can specify both an array of variables and values, in which case the values will be assigned to the corresponding index in the variables array (e.g. the first entry of the value array will be assigned to the first entry of the variable array). In this scenario the user can also choose to evaluate the gradient only w.r.t. certain variables. In any case, the number of values must match the number of specified variables, or if no variables are specified, to the overall amount of variables. 

```jldoctest; setup =:(model = Model(); @variable(model,x); @variable(model,y); @objective(model, Max, x*y - x^2 - y^2))
julia> obj_der(model,[2])
ERROR: There must be as many variables as values
[...]

julia> obj_der(model,[x,y],[1])
ERROR: There must be as many variables as values
[...]
```
The error in the first case results from the fact that the there are two variables registered in the model, while the user tries to only supply one value. The opposite holds true for the second scenario: Here, the user defines two variables, but only supplies one value. 

The `obj_der!` function extends the previously described functionalities by registering the retrieved results within the JuMP model itself. 
```jldoctest example2
julia> model = Model();

julia> @variable(model,x);

julia> @variable(model,y);

julia> @objective(model, Max, x*y - x - y);

julia> obj_der!(model)
2-element Array{GenericAffExpr{Float64,VariableRef},1}:
 y - 1
 x - 1

julia> model[:obj_der]
2-element Array{GenericAffExpr{Float64,VariableRef},1}:
 y - 1
 x - 1

julia> obj_der!(model,[x],:examplename)
1-element Array{GenericAffExpr{Float64,VariableRef},1}:
 y - 1

julia> model[:examplename]
1-element Array{GenericAffExpr{Float64,VariableRef},1}:
 y - 1
```
The default name for the registered gradient is `obj_der`. The result can then be called by querying `model[:obj_der]` (where model is the name of the JuMP model). By specifying a different name in the function arguments we can register the result under a different name, and make it callable under that name. This was also necessary in the above example, since `obj_der!` would have tried to register the default name again, which would have resulted in an error: 
```jldoctest example2
julia> obj_der!(model,[x])
ERROR: You are trying to register the name obj_der in the model, which is already occupied. If you have not specified the name obj_der, then the function you have called assumes obj_der to be the default name. To overcome this issue specify a different name in the function arguments or unregister the name obj_der from the model.
[...]
```
## Constraints 
### Overview 
Different dispatches of the function `con_der` are available to compute the gradient of a specified constraint registered in the JuMP-model `model`. The basic functionalities are the following: 
```@docs 
con_der
```
Results retrieved by these methods can also be registered in the model by calling the function `obj_der!`. 
```@docs
con_der!
```


### Examples
The function `con_der` offers the same functionality as `obj_der`, which shall be be showcased in the following. Consider a constraint which imposes any solution to lie within the unit circle. For the sake of the example we will give this constraint the reference `circle`.
```jldoctest example3
julia> model = Model();

julia> @variable(model,x);

julia> @variable(model,y);

julia> @constraint(model, circle, x^2 + y^2 <= 1);

julia> con_der(model,circle)
2-element Array{GenericAffExpr{Float64,VariableRef},1}:
 2 x
 2 y

julia> con_der(model,circle,[x])
1-element Array{GenericAffExpr{Float64,VariableRef},1}:
 2 x

julia> con_der(model,circle,[3,2])
2-element Array{GenericAffExpr{Float64,VariableRef},1}:
 6
 4

julia> con_der(model,circle,[y],[2])
1-element Array{GenericAffExpr{Float64,VariableRef},1}:
 4
```
Note again that the above function calls of `con_der` are independant of each other. Furthermore, the same remarks hold as for `obj_der`, namely that the values for variables are evaluated according to the order in which the variables are loaded in the model (unless specified otherwise by the user), and that the number of variable values must match the number of (specified) variables. Otherwise, an error will be thrown: 
```jldoctest example3
julia> con_der(model,circle,[1])
ERROR: There must be as many variables as values
[...]

julia> con_der(model,circle,[x,y],[2])
ERROR: There must be as many variables as values
[...]
```
The error in the first case results from the fact that the there are two variables registered in the model, while the user tries to only supply one value. The opposite holds true for the second scenario: Here, the user defines two variables, but only supplies one value. 

The `con_der!` function extends the previously described functionalities by registering the retrieved results within the JuMP model itself. 
```jldoctest example3
julia> con_der!(model, circle)
2-element Array{GenericAffExpr{Float64,VariableRef},1}:
 2 x
 2 y
 
julia> model[:con_der]
2-element Array{GenericAffExpr{Float64,VariableRef},1}:
 2 x
 2 y

julia> con_der!(model,circle,:circle_der)
2-element Array{GenericAffExpr{Float64,VariableRef},1}:
 2 x
 2 y

julia> model[:circle_der]
2-element Array{GenericAffExpr{Float64,VariableRef},1}:
 2 x
 2 y

julia> con_der!(model,circle,[x], :circle_der_x)
1-element Array{GenericAffExpr{Float64,VariableRef},1}:
 2 x

julia> model[:circle_der_x]
1-element Array{GenericAffExpr{Float64,VariableRef},1}:
 2 x
```
 The default name for the registered gradient is `con_der`. The result can then be called by querying `model[:con_der]` (where model is the name of the JuMP model). By specifying a different name in the function arguments we can register the result under a different name, and make it callable under that name. Note how in the above example we registered the results with three different names: First with the default name `con_der`, then with the name `circle_der` and lastly with `circle_der_x`. If we would have attempted to register the default name twice, an error would have been thrown: 
```jldoctest example3
julia> con_der!(model,circle,[x])
ERROR: You are trying to register the name con_der in the model, which is already occupied. If you have not specified the name con_der, then the function you have called assumes con_der to be the default name. To overcome this issue specify a different name in the function arguments or unregister the name con_der from the model.
[...]
```
A similar message would have been thrown if it were attempted to register the same user specified name twice: 
```jldoctest example3
julia> con_der!(model,circle,[y],:circle_der)
ERROR: You are trying to register a name which is already occupied in the model. Choose a different one or unregister the existing.
```
As the error messages suggests it is not strictly necessary to choose a different name for every result registration, one could also unregister any pre-existig name: 
```jldoctest example3
julia> unregister(model,:circle_der)

julia> con_der!(model,circle,[y],:circle_der)
1-element Array{GenericAffExpr{Float64,VariableRef},1}:
 2 y
```
