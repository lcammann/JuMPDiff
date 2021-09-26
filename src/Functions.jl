### Jacobians


"""
model_jac(model)
This function takes as input a JuMP model and returns the Jacobian matrix of the model constraints w.r.t. all variables within the model. 
"""
function model_jac(model::Model)
vars = all_variables(model)
model_jac(model,vars)
end 

""" 
model_jac(model,vars)
This function takes as input a JuMP model and an array of associated variables and returns the Jacobian matrix of the model constraints 
w.r.t. the passed variables. 
"""
function model_jac(model::Model, vars::Array{JuMP.VariableRef})
    #Access all constraints and variables  
    ctypes, cfuncs, crhs, ncto,nct = get_constraint_functions(model)
    allvars = all_variables(model)
    #Pre-allocate output array
    nv = length(vars);          #Number of variables
    jac = Array{GenericAffExpr{Float64,VariableRef}}(undef,(ncto+1,nv))
    #Load partial derivatives
    local ncto;  
    ncto = 1;
    jac[1,:] = reshape(obj_der(model,vars),1,nv) #Load derivative of objective as first entry - Also checks NL expressions!
    for k in 1:nct
        #How many constraints per constraint type ? 
        tmp = MOI.get(model, MOI.ListOfConstraintIndices{ctypes[k][1],ctypes[k][2]}());
        nc = length(tmp);
          for j in 1:nc #Loop over different constraints, per type
            ncto += 1;  #Counter for total number of constraints, independant of type
               for i in 1:nv #Loop of variables, per constraints, per type
                jac[ncto,i] = differentiate_function_new(cfuncs[k][j], vars[i],allvars)   
               end 
           end
   end

    return jac
end 

""" 
model_jac!(model)
This function takes as input a JuMP model and computes its Jacobian. The Jacobian is registered as expression within the model with the default name 
"Jacobian".
"""
function model_jac!(model::Model, name = :Jacobian)
vars = all_variables(model)
model_jac!(model,vars,name)
end 

"""
model_jac!(model,vars)
This function takes as input a JuMP model and an array of associated variables and computes the Jacobian w.r.t. the specified variables. The Jacobian
is registered as expression within the model with the default name "Jacobian". 
"""
function model_jac!(model::Model, vars::Array{JuMP.VariableRef}, name = :Jacobian)
checkregistration(model,name)
jac = model_jac(model,vars)
model.obj_dict[name] = jac;
end 


"""
model_jac(model,vals)
This function takes as input a JuMP model and an array of associated values and returns the Jacobian of the model, evaluated at the specified point.
The number of specified values must match exactly the number of variables within the model. 
"""
function model_jac(model::Model, vals::Array{<:Number})
#Check inputs
vars = all_variables(model)
model_jac(model,vars,vals)
end 


""" 
model_jac(model,vars,vals)
This function takes as input a JuMP model and an array of associated variables and returns the Jacobian matrix of the model constraints 
w.r.t. the passed variables, evaluated at the point specified by vals. 
"""
function model_jac(model::Model, vars::Array{JuMP.VariableRef}, vals::Array{<:Number})
    #Check inputs 
    checklength(vars, vals)
    #Access all constraints and variables  
    ctypes, cfuncs, crhs, ncto,nct = get_constraint_functions(model)
    allvars = all_variables(model)
    #Pre-allocate arrays
    nv = length(vars);          #Number of variables
    jac  = Array{GenericAffExpr{Float64,VariableRef}}(undef,(ncto+1,nv));
    vals = convert(Array{Float64},vals);
    #Create allvals array 
    allvals = convert(Array{GenericAffExpr{Float64,VariableRef}}, allvars);
    for i = 1:nv 
    pos   = findall(isequal(vars[i]), allvars);
    allvals[pos] .= vals[i];
    end 
    #Load partial derivatives
    local ncto;  
    ncto = 1;
    jac[1,:] = reshape(obj_der(model,vars,vals),1,nv) #Load derivative of objective as first entry - Also checks NL expressions!
    for k in 1:nct
        #How many constraints per constraint type ? 
        tmp = MOI.get(model, MOI.ListOfConstraintIndices{ctypes[k][1],ctypes[k][2]}());
        nc = length(tmp);
          for j in 1:nc #Loop over different constraints, per type
            ncto += 1;  #Counter for total number of constraints, independant of type
               for i in 1:nv #Loop of variables, per constraints, per type
                jac[ncto,i] = differentiate_function_new(cfuncs[k][j],vars[i],allvars,vals[i],allvals)   
               end 
           end
   end
 return jac
end 

### Hessians


"""
model_hess(model)
This function returns the Hessian matrix of the objective function specified in model w.r.t to all variables registered in that model. 
"""
function model_hess(model::Model)
vars = all_variables(model)  
model_hess(model::Model,vars::Array{JuMP.VariableRef})
end 

"""
model_hess(model, vars)
This function returns the Hessian matrix of the objective function specified in model w.r.t to the variables specified by vars. 
"""
function model_hess(model::Model,vars::Array{JuMP.VariableRef})
#Check for non-linear expressions 
checknlp(model)
#Load/pre-allocate 
allvars  = all_variables(model); 
obj  = get_objective_function(model); #Note:Currently converts VariableRef objective to scalar affine
nv   = length(vars);
grad = Array{Any,1}(undef,nv);
hess = Array{Float64}(undef,(nv,nv));
#Hessian empty if objective is single variable 
if objective_function_type(model) == VariableRef
    hess = zeros(nv,nv);
    return hess 
end 
#Compute the gradient vector
for i in 1:nv
grad[i] =  differentiate_function_new(obj,vars[i],allvars);
#Assemble Hessian
for j in 1:nv 
    hess[i,j] = 0;
    for k in 1:length(grad[i].terms.keys) #Loop over all variables contained in gradient entry of 
        if grad[i].terms.keys[k] == vars[j]
        hess[i,j] += grad[i].terms.vals[k,1];  
        end 
    end 
end 
end 
return hess
end 

"""
model_hess!(model;name)
This function returns the Hessian matrix of the objective function specified in the model w.r.t. all variables in that model, and additionally registers the results in the model itself.
"""
function model_hess!(model::Model, name = :Hessian)
vars = all_variables(model);
model_hess!(model,vars,name)
end 

"""
model_hess!(model,vars;name)
This function returns the Hessian matrix of the objective fnction specified in model w.r.t. the variables specified by vars, and additionally registers the results in the model itself.
"""
function model_hess!(model::Model,vars::Array{JuMP.VariableRef}, name = :Hessian)
checkregistration(model,name)
hess = model_hess(model,vars);
model.obj_dict[name] = hess
end 

### Constraints


"""
con_der(model,con)
This function takes as input a JuMP model and a constraint reference and returns the gradient vector of the constraint w.r.t. all variables registered in 
the model. 
"""
function con_der(model::Model, con::ConstraintRef)
#Load all variables
vars  = all_variables(model);
#Return derivative
con_der(model,con,vars)
end

"""
con_der(model,con,vars)
This function takes as input a JuMP model, a constraint reference, and an array of variable references and returns the gradient vector of the constraint
w.r.t. all variables specified in vars. 
"""
function con_der(model::Model, con::ConstraintRef, vars::Array{JuMP.VariableRef})
#Check non-linear expressions 
checknlp(model)
#Load constraint function and aux. variables
allvars = all_variables(model);
nv = length(vars);
der = Array{GenericAffExpr{Float64,VariableRef},1}(undef,nv);
bkend = backend(model);
con_func = MOI.get(bkend, MOI.ConstraintFunction(), con.index);
#Return derivative
for i in 1:nv
    der[i] = differentiate_function_new(con_func, vars[i], allvars)
    end 
    return der
end 

"""
con_der(model,con,vals)
This function takes as input a JuMP model, a constraint reference and an array of values and evaluates the derivative of the constraint specified at the values in vals.
"""
function con_der(model::Model, con::ConstraintRef, vals::Array{<:Number})
#Load variables and check input lengths
vars = all_variables(model);
con_der(model,con,vars,vals)
end  


"""
con_der(model,con,vars,vals)
This function takes as input a JuMP model, a constraint reference, an array of variable references and an array of values and evaluates the derivative w.r.t. the specified variables at the specified values.
"""
function con_der(model::Model, con::ConstraintRef, vars::Array{JuMP.VariableRef}, vals::Array{<:Number})
#Check input lengths
checklength(vars,vals)
#Check non-linear expressions
checknlp(model)
#Load relevant variables 
vals = convert(Array{Float64},vals);
allvars = all_variables(model);
nv = length(vars);
der = Array{GenericAffExpr{Float64,VariableRef},1}(undef,nv);
bkend = backend(model);
con_func = MOI.get(bkend, MOI.ConstraintFunction(), con.index);
#Create array "allvals" containing user specified values and all variables
allvals = convert(Array{GenericAffExpr{Float64,VariableRef}}, allvars);
for i = 1:nv 
    pos   = findall(isequal(vars[i]), allvars);
    allvals[pos] .= vals[i];
end 
#Compute derivative 
for i in 1:nv
    der[i] = differentiate_function_new(con_func,vars[i],allvars,vals[i],allvals);
    end 
return der 
end 

"""
con_der!(model,con;name)
This function takes as input a JuMP model, a constraint reference, and optionally a name. It registers the vector of partial derivatives of the constraint within the model. If the user does not specify a value for name, the default name is con_der. Otherwise it is the user specified name.
"""
function con_der!(model::Model, con::ConstraintRef, name = :con_der)
#Load all variables 
vars = all_variables(model);
#Return derivative
con_der!(model, con, vars, name)
end 

"""
con_der!(model,con,vars;name)
This function takes as input a JuMP model, a constraint reference, an array of variables and optionally a symbol. It registers the derivate of the specified constraint at the specified variables within the model. If the user does not specify a value for name, the default name is con_der. Otherwise it is the user specified name.
"""
function con_der!(model::Model, con::ConstraintRef,vars::Array{JuMP.VariableRef}, name = :con_der)
checkregistration(model,name)
der = con_der(model,con,vars)
model.obj_dict[name] = der;
end 

###Objective 


"""
obj_der(model)
This function takes as input a JuMP model and returns the gradient vector of the constraint w.r.t. all variables registered in 
the model. 
"""
function obj_der(model::Model)
#Step 1: Load objective and variables
vars  = all_variables(model);
#Step 2: Return derivative
obj_der(model,vars)
end

"""
obj_der(model,vars)
This function takes as input a JuMP model and an array of variable references and returns the gradient vector of the objective function
w.r.t. all variables specified in vars. 
"""
function obj_der(model::Model, vars::Array{JuMP.VariableRef})
#Step 1: Check non-linear expressions
checknlp(model)
#Step 2: Load objective function and aux. variables
obj     = get_objective_function(model);
allvars = all_variables(model);
nv = length(vars);
der = Array{GenericAffExpr{Float64,VariableRef},1}(undef,nv);
#Step 3: Return derivative
for i in 1:nv
    der[i] = differentiate_function_new(obj, vars[i], allvars)
    end 
    return der
end 

"""
obj_der(model,vals)
This function takes as input a JuMP model and an array of values and evaluates the derivative of the objective function specified at the values in vals.
"""
function obj_der(model::Model, vals::Array{<:Number})
#Load variables and check input lengths
vars = all_variables(model);
obj_der(model,vars,vals)
end  


"""
obj_der(model,vars,vals)
This function takes as input a JuMP model, an array of variable references and an array of values and evaluates the derivative of the objective function w.r.t. the specified variables at the specified values.
"""
function obj_der(model::Model,  vars::Array{JuMP.VariableRef}, vals::Array{<:Number})
#Check input lengths
checklength(vars,vals)
#Check non-linear expressions 
checknlp(model)
#Load objective and  relevant variables 
obj     = get_objective_function(model);
vals    = convert(Array{Float64},vals)
allvars = all_variables(model);
nv = length(vars);
der = Array{GenericAffExpr{Float64,VariableRef},1}(undef,nv);
#Create array "allvals" containing user specified values and all variables
allvals = convert(Array{GenericAffExpr{Float64,VariableRef}}, allvars);
for i = 1:nv 
    pos   = findall(isequal(vars[i]), allvars);
    allvals[pos] .= vals[i];
end 
#Compute derivative 
for i in 1:nv
    der[i] = differentiate_function_new(obj,vars[i],allvars,vals[i],allvals); 
    end 
return der 
end 

"""
obj_der!(model;name)
This function takes as input a JuMP model and optionally a name. It registers the gradient of the objective function within the model under the specified name. If the user does not specify a name, the default name is obj_der. 
"""
function obj_der!(model::Model, name = :obj_der)
#Step 1: Load all variables 
vars = all_variables(model);
#Step 2: Return derivative
obj_der!(model, vars, name)
end 

"""
obj_der!(model,vars;name)
This function takes as input a JuMP model, an array of variable references and optionally a name. It registers the gradient of the objective function w.r.t the specified variables within the model under the specified name. If the user does not specify a name, the default name is obj_der.
"""
function obj_der!(model::Model, vars::Array{JuMP.VariableRef}, name = :obj_der)
checkregistration(model,name)
der = obj_der(model,vars)
model.obj_dict[name] = der;
end 


### Background functions 

#Getting constraints from model
function get_constraint_functions(model::Model)
    local ncto
    bkend = backend(model);
    con_types = MOI.get(model, MOI.ListOfConstraints());
    nct = length(con_types);    
    con_func = Array{Array{MOI.AbstractScalarFunction,1},1}(undef, nct);
    con_rhs = Array{Array{MOI.AbstractSet,1},1}(undef, nct);
    ncto = 0;
    for i in 1:nct
        tmp = MOI.get(model, MOI.ListOfConstraintIndices{con_types[i][1],con_types[i][2]}());
        nc = length(tmp);
        con_func[i] = Array{MOI.AbstractScalarFunction,1}(undef, nc);
        con_rhs[i] = Array{MOI.AbstractSet,1}(undef, nc)
        for j in 1:nc
            ncto += 1;
            con_func[i][j] = MOI.get(bkend, MOI.ConstraintFunction(), tmp[j])
            con_rhs[i][j] = MOI.get(bkend, MOI.ConstraintSet(), tmp[j])
        end
    end
    return con_types, con_func, con_rhs, ncto,nct;
end


#Differentiating single terms
function differentiate_term(term::MOI.ScalarAffineTerm{Float64}, var_id::MOI.VariableIndex)
    return term.variable_index == var_id ? term.coefficient : 0.;
end


#New function for differentiating quadratic terms
function differentiate_function_new(func::MOI.ScalarQuadraticFunction{Float64}, var::VariableRef,vars::Array{JuMP.VariableRef})
    local  mults #Bind mults outisde of loop
    nq       = length(func.quadratic_terms);
    constant = 0.;
    mults    = Array{Any}(undef, (1,nq));
    for term in func.affine_terms 
        constant += differentiate_term(term, index(var));
    end
    for i in 1:nq
        coeff    = 0.;
        varl     = var; #Helper variable to not overwrite var itself in the loop - There should be a better way to do this
        term = func.quadratic_terms[i];
        #Simple quadratic
        if term.variable_index_1 == term.variable_index_2 == index(var)
            coeff = term.coefficient;
        #Bilinear - Case 1
        elseif term.variable_index_1 != term.variable_index_2 && term.variable_index_2 == index(var)
            coeff = term.coefficient; 
            pos   = findall(isequal(term.variable_index_1), index.(vars)); #Position in vars vector for unknown variable 
            try 
            varl   = vars[pos[1,1]];
            catch #If user only specifies one of the bilinear variables, varl = 0
            varl   = 0.;
            end 
        #Bilinear - Case 2
        elseif term.variable_index_2 != term.variable_index_1 && term.variable_index_1 == index(var)
            coeff = term.coefficient; 
            pos   = findall(isequal(term.variable_index_2), index.(vars)); #Position in vars vector for unknown variable 
            try
            varl  = vars[pos[1,1]];
            catch #If user only specifies one of the bilinear variables, varl = 0
            varl  = 0.;
            end 
        end 
    mults[i] = coeff*varl
    end 
return sum(mults) + constant
end 

#New function for differentiating quadratic terms with values
function differentiate_function_new(func::MOI.ScalarQuadraticFunction{Float64}, var::VariableRef, vars::Array{JuMP.VariableRef}, val::Float64, vals::Array{GenericAffExpr{Float64,VariableRef}})
    local mults #Bind mults outisde of loop
    nq       = length(func.quadratic_terms);
    constant = 0.;
    mults    = Array{Any}(undef, (1,nq));
    for term in func.affine_terms 
        constant += differentiate_term(term, index(var));
    end
    for i in 1:nq
        coeff    = 0.;
        thisval  = val; #Helper variable to not overwrite val in the loop 
        term = func.quadratic_terms[i];
       #Simple quadratic
       if term.variable_index_1 == term.variable_index_2 == index(var)
        coeff = term.coefficient;
    #Bilinear - Case 1
    elseif term.variable_index_1 != term.variable_index_2 && term.variable_index_2 == index(var)
        coeff = term.coefficient; 
        pos   = findall(isequal(term.variable_index_1), index.(vars)); #Position in vars vector for unknown variable 
        thisval   = vals[pos[1,1]];
    #Bilinear - Case 2
    elseif term.variable_index_2 != term.variable_index_1 && term.variable_index_1 == index(var)
        coeff = term.coefficient; 
        pos   = findall(isequal(term.variable_index_2), index.(vars)); #Position in vars vector for unknown variable 
        thisval   = vals[pos[1,1]];
    end 
mults[i] = coeff*thisval
end 
return sum(mults) + constant
end 


#Takes a linear function and a model variable and returns the derivative w.r.t. the 
#given variable
function differentiate_function_new(func::MOI.ScalarAffineFunction{Float64}, var::VariableRef,vars::Array{JuMP.VariableRef},val=0.0, vals=[0.0])
    constant = 0.;
    for term in func.terms
        constant += differentiate_term(term, index(var));
    end
    return constant;
end

#Takes a function composed of a single variable and a model variable and returns 
#the derivative w.r.t. the given variable
function differentiate_function_new(func::MOI.SingleVariable, var::VariableRef,vars::Array{JuMP.VariableRef},val=0.0, vals=[0.0])
    return func.variable == index(var) ? 1. : 0.;
end

#Checks whether a user supplied Symbol is already registered in the model. If true, returns an error message. 
function checkregistration(model::Model, name::Symbol)
    if haskey(model.obj_dict, name) == true && name == :Jacobian
        error("You are trying to register the name Jacobian in the model, which is already occupied. If you have not specified the name Jacobian, then the function you have called assumes Jacobian to be the default name. To overcome this issue specify a different name in the function arguments or unregister the name Jacobian from the model.")
    elseif haskey(model.obj_dict, name) == true && name == :Hessian
        error("You are trying to register the name Hessian in the model, which is already occupied. If you have not specified the name Hessian, then the function you have called assumes Hessian to be the default name. To overcome this issue specify a different name in the function arguments or unregister the name Hessian from the model.")
    elseif haskey(model.obj_dict, name) == true && name == :con_der 
        error("You are trying to register the name con_der in the model, which is already occupied. If you have not specified the name con_der, then the function you have called assumes con_der to be the default name. To overcome this issue specify a different name in the function arguments or unregister the name con_der from the model.")
    elseif haskey(model.obj_dict, name) == true && name == :obj_der 
        error("You are trying to register the name obj_der in the model, which is already occupied. If you have not specified the name obj_der, then the function you have called assumes obj_der to be the default name. To overcome this issue specify a different name in the function arguments or unregister the name obj_der from the model.")
    elseif haskey(model.obj_dict, name) == true
        error("You are trying to register a name which is already occupied in the model. Choose a different one or unregister the existing.")
end 
end 

#Checks whether length of variable and value vectors matches 
function checklength(vars::Array{JuMP.VariableRef}, vals::Array{<:Number})
    if length(vars) != length(vals)
        error("There must be as many variables as values")
    end 
end

#Cheks whether non-linear expressions are present
function checknlp(model::Model)
    if model.nlp_data !== nothing 
        println("You have non-linear expressions registered in your model, which are currently not supported by JuMP-Diff. The function you have called will ignore these expressions and proceed.")
    end 
end
    

#Returns the objective function of the given model as a MathOptInterface object 
function get_objective_function(model::Model)
    if objective_function_type(model) == GenericAffExpr{Float64,VariableRef}
        return MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}());
    elseif objective_function_type(model) == GenericQuadExpr{Float64,VariableRef}
        return MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}}());
    elseif objective_function_type(model) == VariableRef
        return MOI.get(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}()); #Should be sensible
    else
        error("Objective function is not affine nor quadratic.");
    end
end
