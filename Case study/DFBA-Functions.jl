function DFBA_SEQUENTIAL(z0,X0;kwargs...)
     t   = get(kwargs, :t, Setup[:t]);
     N   = get(kwargs, :N, Setup[:N]);
     St  = get(kwargs, :St, Setup[:St]);
     kla = get(kwargs, :kla, Setup[:kla]);
     Km  = get(kwargs, :Km, Setup[:Km]);
     ΔT  = t[end]/N;
          
     zval = Array{Float64}(undef,(3,N));
     vval = Array{Float64}(undef,(4,N));
     Xval = Array{Float64}(undef,(1,N));
          
     zval[:,1] = z0';
     Xval[1]   = X0;
     w         = [1,1,1,1];               #Weights
     
     
     #Loop over optimization with instantaneous objective
     for i in 1:N-1
          model = Model(Ipopt.Optimizer)
          set_optimizer_attributes(model,"print_level"=> 0)
          @variable(model, v[1:4] >= 0);
          @variable(model, z[1:3,1:2] >= 0);
          @variable(model, X[1:2] >= 0);
     
          @constraint(model,z[:,1] .== zval[:,i]);
          @constraint(model,X[1] == Xval[i]);
          @constraint(model, -St[3,:]'*v <= 15);
          @constraint(model, -(St[1,2]*v[2] + St[1,3]*v[3]  + St[1,4]*v[4])*(Km + z[1,2]) <= 10*z[1,2]);
          @constraint(model, z[1,2] == z[1,1] + St[1,:]'*v*ΔT*X[2]);
          @constraint(model, z[2,2] == z[2,1] + St[2,:]'*v*ΔT*X[2]);
          @constraint(model, z[3,2] == z[3,1] + St[3,:]'*v*ΔT*X[2] + kla*(0.21 - z[3,2])*ΔT);
          @constraint(model, X[2]   .== X[1] + sum(v)*X[2]*ΔT);
     
          @objective(model,Max,w'*v);
          optimize!(model)
          println(termination_status(model))
     
          z_val = value.(z);
          zval[:,i+1] = z_val[:,2];
          vval[:,i] = value.(v);
          X_val = value.(X);
          Xval[i+1] = X_val[2];
     end
     
     return vval,zval,Xval
     end 

function DFBA_KKT(z0,X0;kwargs...)
t   = get(kwargs, :t, Setup[:t]);
N   = get(kwargs, :N, Setup[:N]);
St  = get(kwargs, :St, Setup[:St]);
kla = get(kwargs, :kla, Setup[:kla]);
Km  = get(kwargs, :Km, Setup[:Km]);
ΔT  = t[end]/N;

zval = Array{Float64}(undef,(3,N));
vval = Array{Float64}(undef,(4,N));
Xval = Array{Float64}(undef,(1,N));

zval[:,1] = z0';
Xval[1]   = X0;
w         = [1,1,1,1];               #Weights

nvar = 34;          #Number of variables
ncon = 22;          #Number of constraints
neq  = 8;           #Number of equality constraints
nieq = ncon - neq;  #Number of inequality constraints
L   = Array{Any}(undef,ncon+1,1)   	          
der = Array{GenericQuadExpr{Float64,VariableRef}}(undef,nvar,ncon+1)
for i in 1:N-1
model = Model(Ipopt.Optimizer);
set_optimizer_attributes(model,"print_level"=> 0);
@variable(model,v[1:4] >= 0);
@variable(model,z[1:3,1:2] >= 0);
@variable(model, X[1:2] >= 0);
@variable(model,λ[1:8]);
@variable(model,μ[1:14] >= 0);
varlist = all_variables(model);

L[1]  = w'v;

## Add equality constraints
for j in 1:3
L[j+1] = (zval[j,i] - z[j,1]);   #Metabolite concentrations at t
end
L[5] = (Xval[i] - X[1]);          #Biomass concentration at t
L[6] = (z[1,2] - z[1,1] - St[1,:]'*v*ΔT*X[2]);
L[7] = (z[2,2] - z[2,1] - St[2,:]'*v*ΔT*X[2]);
L[8] = (z[3,2] - z[3,1] - St[3,:]'*v*ΔT*X[2] - kla*(0.21 - z[3,2])*ΔT);
L[9] = (X[2]   - X[1] - sum(v)*X[2]*ΔT);

@objective(model,Max,L[1])
der[:,1] = obj_der(model)
global ∇L = der[:,1];
for j in 2:neq+1
@objective(model,Max,L[j])
der[:,j] = λ[j-1]*obj_der(model) #Hacky workaround, for otherwise we are dealing with trilinear terms
global ∇L += der[:,j]
end

#Add partials w.r.t. multipliers (i.e., the constraints)
for j in 1:length(λ)
global ∇L[j+nvar-length(μ)-length(λ)] = L[j+1]
end


## Add inequality constraints
L[10] = (St[1,2]*v[2] + St[1,3]*v[3]  + St[1,4]*v[4])*(Km + z[1,2]) + 10*z[1,2];
L[11] = St[3,:]'*v + 15;

# Inequality constraints regarding single variables
for j in 1:nvar-length(μ)-length(λ)
     L[nvar-length(μ)-length(λ)-1+j] = varlist[j]
end

for j in 1:length(μ)
     @objective(model,Max,L[length(λ)+1+j])
     der[:,length(λ)+1+j] = μ[j]*obj_der(model)
     global ∇L +=der[:,length(λ)+1+j]
end

#Add constraints to the model
@constraint(model, EQC[i in 1:nvar-length(μ)], ∇L[i] == 0);
@constraint(model, IEQC[i in neq+2:ncon], L[i] >= 0);


@NLobjective(model,Min, sum((μ[i]*L[neq+1+i]) for i in 1:nieq))
optimize!(model)

#Save values for next iteration
z_val = value.(z);
zval[:,i+1] = z_val[:,2];
vval[:,i] = value.(v);
X_val = value.(X);
Xval[i+1] = X_val[2];
end

l  = @layout[a;b;c]
p1 = plot(LinRange(0,t[end],N), value.(vval)', title = "Sequential, KKT", label = ["v1" "v2" "v3" "v4"], xlabel = "t [h]", ylabel = "[mmol/(gdw*hr)]")
p2 = plot(LinRange(0,t[end],N), value.(zval)', label = ["Glc" "Ac" "O2"], xlabel = "t [h]", ylabel = "[mmol]" )
p3 = plot(LinRange(0,t[end],N), value.(Xval)', label = "Biom.", xlabel = "t [h]", ylabel = "[g/l]")

return vval,zval,Xval
end

function DFBA_SIMULTANEOUS(z0,X0;kwargs...)
t   = get(kwargs, :t, Setup[:t]);
N   = get(kwargs, :N, Setup[:N]);
St  = get(kwargs, :St, Setup[:St]);
kla = get(kwargs, :kla, Setup[:kla]);
Km  = get(kwargs, :Km, Setup[:Km]);
ncp = get(kwargs, :ncp, Setup[:ncp]);
ΔT  = t[end]/N;

model = Model(Ipopt.Optimizer)
@variable(model,v[1:N,1:4] >= 0);
@variable(model,z[1:N,1:3,1:ncp+1] >= 0);
@variable(model,X[1:N, 1:ncp+1] >= 0);
@variable(model, 0 <= D[1:N] <= 1);
zin = z0;
zin[3] = 0;

#Equality constraints
@constraint(model,ode_g[i in 1:N,k in 2:ncp+1],sum(z[i,1,l]*adot[l,k] for l in 1:ncp+1) - ΔT*(sum(St[1,r]*v[i,r] for r in 1:4)*X[i,k] + D[i]*(zin[1] - z[i,1,k])) == 0);
@constraint(model,ode_a[i in 1:N,k in 2:ncp+1],sum(z[i,2,l]*adot[l,k] for l in 1:ncp+1) - ΔT*(sum(St[2,r]*v[i,r] for r in 1:4)*X[i,k] + D[i]*(zin[2] - z[i,2,k])) == 0);
@constraint(model,ode_o[i in 1:N,k in 2:ncp+1],sum(z[i,3,l]*adot[l,k] for l in 1:ncp+1) - ΔT*(sum(St[3,r]*v[i,r] for r in 1:4)*X[i,k] + kla*(0.21 - z[i,3,k]) + D[i]*(zin[3] - z[i,3,k])) == 0);
@constraint(model,ode_x[i in 1:N,k in 2:ncp+1],sum(X[i,l]*adot[l,k] for l in 1:ncp+1) - ΔT*(sum(v[i,r] for r in 1:4)*X[i,k] + D[i]*(-X[i,k])) == 0);

#Continuity constraints
@constraint(model,cont_z[i in 2:N, j in 1:3], z[i,j,1] == z[i-1,j,ncp+1]);
@constraint(model,cont_x[i in 2:N], X[i,1] == X[i-1,ncp+1]);

#Inequality constraints
@constraint(model,maxO2[i in 1:N], St[3,:]'*v[i,:] + 15 >= 0);
@constraint(model,MM[i in 1:N], (St[1,2]*v[i,2] + St[1,3]*v[i,3]  + St[1,4]*v[i,4])*(Km + z[i,1,ncp+1]) + 10*z[i,1,ncp+1] >= 0);


#Objective
@objective(model,Max,sum(X[i,r]/X0 for i in 1:N for r in 1:4)) #Equivalent to maximizing fluxes

#Initial values
fix.(z[1,:,1], z0, force = true);
fix(X[1,1], X0, force = true);
fix.(D[1:N], 0, force = true);
optimize!(model)

vval = value.(v)
zval = value.(z[:,:,1])
Xval = value.(X[:,1])

return vval,zval,Xval
end

function DFBA_COMPARISON(z0,X0)
vval,zval,Xval = DFBA_KKT(z0, X0 ; N = 100, t = 10)
p1 = plot(LinRange(0,10,100), value.(vval)', title = "Sequential, KKT", label = ["v1" "v2" "v3" "v4"], xlabel = "t [h]", ylabel = "[mmol/(gdw*hr)]")
p2 = plot(LinRange(0,10,100), value.(zval)', label = ["Glc" "Ac" "O2"], xlabel = "t [h]", ylabel = "[mmol]" )
p3 = plot(LinRange(0,10,100), value.(Xval)', label = "Biom.", xlabel = "t [h]", ylabel = "[g/l]")
     
vval,zval,Xval = DFBA_SEQUENTIAL(z0, X0 ; N = 100, t=10)
p4 = plot(LinRange(0,10,100), value.(vval)', title = "Sequential",label = ["v1" "v2" "v3" "v4"], xlabel = "t [h]", ylabel = "[mmol/(gdw*hr)]")
p5 = plot(LinRange(0,10,100), value.(zval)', label = ["Glc" "Ac" "O2"], xlabel = "t [h]", ylabel = "[mmol]" )
p6 = plot(LinRange(0,10,100), value.(Xval)', label = "Biom.", xlabel = "t [h]", ylabel = "[g/l]")
     
vval,zval,Xval = DFBA_SIMULTANEOUS(z0, X0 ; N = 100, t=10)
p7 = plot(LinRange(0,10,100), value.(vval), title = "Simultaneously",label = ["v1" "v2" "v3" "v4"], xlabel = "t [h]", ylabel = "[mmol/(gdw*hr)]")
p8 = plot(LinRange(0,10,100), value.(zval), label = ["Glc" "Ac" "O2"], xlabel = "t [h]", ylabel = "[mmol]" )
p9 = plot(LinRange(0,10,100), value.(Xval), label = "Biom.", xlabel = "t [h]", ylabel = "[g/l]")
     
DFBAplot = plot(p1,p4,p7,p2,p5,p8,p3,p6,p9, layout = (3,3))
return DFBAplot
end 