using JuMPDiff
model = Model();
@variable(model,x);
@variable(model,y)
@constraint(model, x - y == 5);
@constraint(model, y^2 + x^2 >= 0);
@constraint(model,3*x^2 >= 5*x); 
@constraint(model,x*y == 10);    #Constraints of this type don't work yet
model_pdv(model,[x,y])
