" " " 
model_pdv(model,vars)

This function takes as input a JuMP model and a set of associated variables and returns the partial 
derivatives of the model constraints w.r.t. to the passed variables. 
" " "

function model_pdv(model::Model, vars::Array{JuMP.VariableRef})
    # Access the constraints 
    ctypes, cfuncs, crhs, ncto,nct = get_constraint_functions(model)
    # Pre-allocate output array
    nv = length(vars);          #Number of variables
    pdv = Array{Any}(undef,(ncto,nv))
    # Load partial derivatives
    local ncto;  
    ncto = 0;
    for k in 1:nct
        # How many constraints per constraint type ? 
        tmp = MOI.get(model, MOI.ListOfConstraintIndices{ctypes[k][1],ctypes[k][2]}());
        nc = length(tmp);
          for j in 1:nc # Loop over different constraints, per type
            ncto += 1;  # Counter for total number of constraints, independant of type
               for i in 1:nv # Loop of variables, per constraints, per type
                pdv[ncto,i] = differentiate_function_new(cfuncs[k][j], vars[i])   
               end 
           end
   end

    return pdv
    # Print results
    print("ctypes =", ctypes)
    print("cfuncs =",  cfuncs)
    print("crhs =", crhs)
end 

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
# returns the derivative of an affine term
function differentiate_term(term::MOI.ScalarAffineTerm{Float64}, var_id::MOI.VariableIndex)
    return term.variable_index == var_id ? term.coefficient : 0.;
end

# returns the derivative of a quadratic term
function differentiate_term(term::MOI.ScalarQuadraticTerm{Float64}, var_id::MOI.VariableIndex)
    if term.variable_index_1 == var_id
        return MOI.ScalarAffineTerm(term.variable_index_2 == 
                var_id ? 2.0*term.coefficient : term.coefficient, term.variable_index_2);
    elseif term.variable_index_2 == var_id
        return MOI.ScalarAffineTerm(term.coefficient, term.variable_index_1);
    else
        return 0;
    end
end
###








### Differentiating functions

##New idea
function differentiate_function_new(func::MOI.ScalarQuadraticFunction{Float64}, var::VariableRef)
    local coeff #Bind coeff outisde of loop
    constant = 0.;
    coeff    = 0.;
    for term in func.affine_terms 
        constant += differentiate_term(term, index(var));
    end
    for term in func.quadratic_terms
        if term.variable_index_1 == index(var)
            coeff = term.coefficient;
        end 
    end 
return coeff*var + constant
end 

# takes a linear function and a model variable and returns the derivative w.r.t. the 
# given variable
function differentiate_function_new(func::MOI.ScalarAffineFunction{Float64}, var::VariableRef)
    constant = 0.;
    for term in func.terms
        constant += differentiate_term(term, index(var));
    end
    return constant;
end

# takes a function composed of a single variable and a model variable and returns 
# the derivative w.r.t. the given variable
function differentiate_function_new(func::MOI.SingleVariable, var::VariableRef)
    return func.variable == index(var) ? 1. : 0.;
end





# takes a quadratic function and a model variable and returns the derivative w.r.t. the 
# given variable
function differentiate_function(func::MOI.ScalarQuadraticFunction{Float64}, var::VariableRef)
    constant = 0.;
    for term in func.affine_terms 
        constant += differentiate_term(term, index(var));
    end
    terms = Array{MOI.ScalarAffineTerm{Float64},1}();
    for term in func.quadratic_terms
        t = differentiate_term(term, index(var))
        if t != 0
            append!(terms, [t]);
        end
    end
    if !isempty(terms)
        return MOI.ScalarAffineFunction(terms, constant);
    else
        return constant; 
    end

end

# takes a linear function and a model variable and returns the derivative w.r.t. the 
# given variable
function differentiate_function(func::MOI.ScalarAffineFunction{Float64}, var::VariableRef)
    constant = 0.;
    for term in func.terms
        constant += differentiate_term(term, index(var));
    end
    return constant;
end

# takes a function composed of a single variable and a model variable and returns 
# the derivative w.r.t. the given variable
function differentiate_function(func::MOI.SingleVariable, var::VariableRef)
    return func.variable == index(var) ? 1. : 0.;
end