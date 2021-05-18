using JuMPDiff
model = Model();
@variable(model,x);
@variable(model,y);
@variable(model,z);
@constraint(model, x - y == 5);
@constraint(model, 3*x^2 + 2*y >= 5*z);
@constraint(model,x*y == 10); 
@objective(model, Max, x^2 + y^2 +z*x);

model_jac(model)
model_jac(model,[x,y])
model_jac(model,[1,1])
model_jac(model,[x,y],[1,1])
model_jac!(model) #Expression registered, but can't easily be called

model_hess(model)
model_hess(model,[x,z]) #Needs reworking 
