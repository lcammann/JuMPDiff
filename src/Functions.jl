"""
model_jac(model)

This function takes an input a JuMP model and returns the Jacobian matrix of the model constraints w.r.t. all variables within the model. 
"""
function model_jac(model::Model)
vars = all_variables(model)
model_jac(model,vars)
end 

""" 
model_jac(model,vars)

This function takes as input a JuMP model and a set of associated variables and returns the Jacobian matrix of the model constraints 
w.r.t. to the passed variables. 
"""
function model_jac(model::Model, vars::Array{JuMP.VariableRef})
    #Access all constraints and variables  
    ctypes, cfuncs, crhs, ncto,nct = get_constraint_functions(model)
    allvars = all_variables(model)
    #Pre-allocate output array
    nv = length(vars);          #Number of variables
    jac = Array{Any}(undef,(ncto,nv))
    #Load partial derivatives
    local ncto;  
    ncto = 0;
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

This function takes as input a JuMP model and computes its Jacobian. The Jaocbian is registered as expression within the model with the default name 
"Jacobian".
"""
function model_jac!(model::Model)
vars = all_variables(model)
model_jac!(model,vars)
end 

"""
model_jac!(model,vars)

This function takes as input a JuMP model and an array of associated variables and computes the Jacobian w.r.t. the specified variables. The Jacobian
is registered as expression within the model with the default name "Jacobian". 
"""
function model_jac!(model::Model, vars::Array{JuMP.VariableRef})
    #Access all constraints and variables  
    ctypes, cfuncs, crhs, ncto,nct = get_constraint_functions(model)
    allvars = all_variables(model)
    #Pre-allocate output array
    nv = length(vars);          #Number of variables
    jac = Array{Any}(undef,(ncto,nv))
    #Load partial derivatives
    local ncto;  
    ncto = 0;
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
@expression(model,Jacobian,jac)
end 


"""
model_jac(model,vals)

This function takes as input a JuMP model and an array of associated values and returns the Jacobian of the model, evaluated at the specified point.
The number of specified values must match exactly the number of variables within the model. 
"""
function model_jac(model::Model, vals::Array{<:Number})
#Check inputs
vars = all_variables(model)
if length(vals) != length(vars)
error("There must be as many variables as values")
end 
model_jac(model,vars,vals)
end 


""" 
model_jac(model,vars,vals)

This function takes as input a JuMP model and a set of associated variables and returns the Jacobian matrix of the model constraints 
w.r.t. to the passed variables, evaluated at the point specified by vals. 
"""
function model_jac(model::Model, vars::Array{JuMP.VariableRef}, vals::Array{<:Number})
    #Check inputs 
    if length(vars) != length(vals)
        error("There must be as many variables as values")
    end 
    #Access all constraints and variables  
    ctypes, cfuncs, crhs, ncto,nct = get_constraint_functions(model)
    allvars = all_variables(model)
    #Pre-allocate arrays
    nv = length(vars);          #Number of variables
    jac  = Array{Any}(undef,(ncto,nv))
    vals = convert(Array{Float64},vals)
    #Load partial derivatives
    local ncto;  
    ncto = 0;
    for k in 1:nct
        #How many constraints per constraint type ? 
        tmp = MOI.get(model, MOI.ListOfConstraintIndices{ctypes[k][1],ctypes[k][2]}());
        nc = length(tmp);
          for j in 1:nc #Loop over different constraints, per type
            ncto += 1;  #Counter for total number of constraints, independant of type
               for i in 1:nv #Loop of variables, per constraints, per type
                jac[ncto,i] = differentiate_function_new(cfuncs[k][j],vars[i],vars,vals[i],vals)   
               end 
           end
   end
 return jac
end 

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
#Load/pre-allocate 
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
grad[i] =  differentiate_function_new(obj,vars[i],vars);
#Assemble 
for j in 1:nv 
    if grad[i].terms.keys == VariableRef[vars[j]]
       hess[i,j] = grad[i].terms.vals[1,1];
    else 
       hess[i,j] = 0;
    end 
end 
end 
return hess
end 

### - Getting constraints from model
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


### - Differentiating single terms
function differentiate_term(term::MOI.ScalarAffineTerm{Float64}, var_id::MOI.VariableIndex)
    return term.variable_index == var_id ? term.coefficient : 0.;
end

### Differentiating functions

##New function for quadratic terms
function differentiate_function_new(func::MOI.ScalarQuadraticFunction{Float64}, var::VariableRef,vars::Array{JuMP.VariableRef})
    local coeff #Bind coeff outisde of loop
    constant = 0.;
    coeff    = 0.;
    for term in func.affine_terms 
        constant += differentiate_term(term, index(var));
    end
    for term in func.quadratic_terms
        #Simple quadratic
        if term.variable_index_1 == term.variable_index_2 == index(var)
            coeff = term.coefficient;
        #Bilinear - Case 1
        elseif term.variable_index_1 != term.variable_index_2 && term.variable_index_2 == index(var)
            coeff = term.coefficient; 
            pos   = findall(isequal(term.variable_index_1), index.(vars)); #Position in vars vector for unknown variable 
            try 
            var   = vars[pos[1,1]];
            catch #If user only specifies one of the bilinear variables
            var   = 0;
            end 
        #Bilinear - Case 2
        elseif term.variable_index_2 != term.variable_index_1 && term.variable_index_1 == index(var)
            coeff = term.coefficient; 
            pos   = findall(isequal(term.variable_index_2), index.(vars)); #Position in vars vector for unknown variable 
            try
            var   = vars[pos[1,1]];
            catch #If user only specifies one of the bilinear variables 
            var   = 0;
            end 
        end 
    end 
return coeff*var + constant
end 

function differentiate_function_new(func::MOI.ScalarQuadraticFunction{Float64}, var::VariableRef, vars::Array{JuMP.VariableRef}, val::Float64, vals::Array{<:Number})
    local coeff #Bind coeff outisde of loop
    constant = 0.;
    coeff    = 0.;
    for term in func.affine_terms 
        constant += differentiate_term(term, index(var));
    end
    for term in func.quadratic_terms
       #Simple quadratic
       if term.variable_index_1 == term.variable_index_2 == index(var)
        coeff = term.coefficient;
    #Bilinear - Case 1
    elseif term.variable_index_1 != term.variable_index_2 && term.variable_index_2 == index(var)
        coeff = term.coefficient; 
        pos   = findall(isequal(term.variable_index_1), index.(vars)); #Position in vars vector for unknown variable 
        try 
        val   = vals[pos[1,1]];
        catch #If user only specifies one of the bilinear variables
        val   = 0;
        end 
    #Bilinear - Case 2
    elseif term.variable_index_2 != term.variable_index_1 && term.variable_index_1 == index(var)
        coeff = term.coefficient; 
        pos   = findall(isequal(term.variable_index_2), index.(vars)); #Position in vars vector for unknown variable 
        try
        val   = vals[pos[1,1]];
        catch #If user only specifies one of the bilinear variables 
        val   = 0;
        end 
    end 
end 
return coeff*val + constant
end 


# Takes a linear function and a model variable and returns the derivative w.r.t. the 
# given variable
function differentiate_function_new(func::MOI.ScalarAffineFunction{Float64}, var::VariableRef,vars::Array{JuMP.VariableRef},val=0.0, vals=[0.0])
    constant = 0.;
    for term in func.terms
        constant += differentiate_term(term, index(var));
    end
    return constant;
end

# takes a function composed of a single variable and a model variable and returns 
# the derivative w.r.t. the given variable
function differentiate_function_new(func::MOI.SingleVariable, var::VariableRef,vars::Array{JuMP.VariableRef},val=0.0, vals=[0.0])
    return func.variable == index(var) ? 1. : 0.;
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