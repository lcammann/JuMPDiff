function Bilevel(z0,X0,flag;kwargs...)
t   = get(kwargs, :t, Setup[:t]);
N   = get(kwargs, :N, Setup[:N]);
St  = get(kwargs, :St, Setup[:St]);
kla = get(kwargs, :kla, Setup[:kla]);
Km  = get(kwargs, :Km, Setup[:Km]);
ΔT  = t[end]/N;
ncp = 3;

#Set up the problem
ic   = 4;                       #Initial conditions
neq  = 4*(N*ncp + (N-1)) + ic; 
nieq = 24*N;
ncon = neq + nieq; 
zin = z0;
zin[3] = 0;

model = Model(Ipopt.Optimizer)
@variable(model,v[1:N,1:4] >= 0);
@variable(model,z[1:N,1:3,1:ncp+1] >= 0);
@variable(model,X[1:N, 1:ncp+1] >= 0);
@variable(model, 0 <= D[1:N] <= 1);
@variable(model, λ[1:neq]);
@variable(model,μ[1:nieq] >= 0);

varlist = all_variables(model);
nvar = length(varlist);


#Equality constraints
@expression(model,ode_g[i in 1:N,k in 2:ncp+1],sum(z[i,1,l]*adot[l,k] for l in 1:ncp+1) - ΔT*(sum(St[1,r]*v[i,r] for r in 1:4)*X[i,k] + D[i]*(zin[1] - z[i,1,k])))
@expression(model,ode_a[i in 1:N,k in 2:ncp+1],sum(z[i,2,l]*adot[l,k] for l in 1:ncp+1) - ΔT*(sum(St[2,r]*v[i,r] for r in 1:4)*X[i,k] + D[i]*(zin[2] - z[i,2,k])));
@expression(model,ode_o[i in 1:N,k in 2:ncp+1],sum(z[i,3,l]*adot[l,k] for l in 1:ncp+1) - ΔT*(sum(St[3,r]*v[i,r] for r in 1:4)*X[i,k] + kla*(0.21 - z[i,3,k]) + D[i]*(zin[3] - z[i,3,k])));
@expression(model,ode_x[i in 1:N,k in 2:ncp+1],sum(X[i,l]*adot[l,k] for l in 1:ncp+1) - ΔT*(sum(v[i,r] for r in 1:4)*X[i,k] + D[i]*(-X[i,k]))); 
@expression(model,cont_z[i in 2:N, j in 1:3], z[i,j,1] - z[i-1,j,ncp+1]);
@expression(model,cont_x[i in 2:N], X[i,1] - X[i-1,ncp+1]);
#Inequality constraints
@expression(model,maxO2[i in 1:N], St[3,:]'*v[i,:] + 15 ); 
@expression(model,MM[i in 1:N], (St[1,2]*v[i,2] + St[1,3]*v[i,3]  + St[1,4]*v[i,4])*(Km + z[i,1,ncp+1]) + 10*z[i,1,ncp+1]);

L = Array{Any}(undef,ncon+1,1);
der = Array{GenericQuadExpr{Float64,VariableRef}}(undef,nvar,ncon+1)
#Add objective function
L[1] = sum(v[i,r] for i in 1:N, r in 1:4) #model[:obj] #sum(X[i,r]/i for i in 1:N for r in 1:4)#model[:obj] 

#Add intial conditions 
for j in 1:3
    L[j+1] = (z0[j] - z[1,j,1]);   #Metabolite concentrations at t
    end 
L[5] = (X0 - X[1,1]);     #Biomass concentration at t

##Set up equality constraints
global j = ic+2;
#Load ode_g
for i in 1:N
    for k in 2:ncp+1
        L[j] = ode_g[i,k];
        global j += 1;     #counter
    end 
end 
#Load ode_a
for i in 1:N
    for k in 2:ncp+1
        L[j] = ode_a[i,k];
        global j += 1;     #counter
    end 
end 
#Load ode_o
for i in 1:N
    for k in 2:ncp+1
        L[j] = ode_o[i,k];
        global j += 1;     #counter
    end 
end
#Load ode_x 
for i in 1:N
    for k in 2:ncp+1
        L[j] = ode_x[i,k];
        global j += 1;     #counter
    end 
end 
#Load continuties 
for i in 2:N
    for k in 1:3 
        L[j] = cont_z[i,k];
        global j += 1; 
    end 
end 
for i in 2:N
    L[j] = cont_x[i]; 
    global j += 1;
end 

@objective(model,Max,L[1])
der[:,1] = obj_der(model)
global ∇L = der[:,1];
for i in 2:neq+1
    @objective(model,Max,L[i])
    der[:,i] = λ[i-1]*obj_der(model) #Hacky workaround, for otherwise we are dealing with trilinear terms
    global ∇L -= der[:,i]
end 


## Set up inequality constraints
#Load Michaelis-Menten 
for i in 1:N
        L[j] = MM[i];
        global j += 1; 
end 
#Load maxO2
for i in 1:N
    L[j] = maxO2[i]; 
    global j += 1;
end

varsnomu = varlist[1:end-length(μ)-length(λ)]   #No inequality constraints for multipliers
for i in 1:length(varsnomu)
    L[j] = varsnomu[i] #Inequality constraints regarding single variable
    global j += 1;
end
L[j:end] = 1.0 .- D;

for i in 1:length(μ)
    @objective(model,Max,L[length(λ)+i])
    der[:,length(λ)+i] = μ[i]*obj_der(model)
    global ∇L -=der[:,length(λ)+i]
end 




#Set up optimization 
@variable(model,-0.7 <= s[1:nieq] <= 0.7); #slack variables

@constraint(model,GRADL[i = 1:nvar-length(μ)-length(λ)], ∇L[i] == 0);
@constraint(model,EQC[i in 1:neq], L[i+1] == 0);
@NLconstraint(model,slack[i in 1:nieq], μ[i]*L[neq+1+i] + s[i]== 0 );
@constraint(model, IEQC[i in neq+2:length(L)], L[i] >= 0);

if flag == 1
@objective(model,Max,sum(D[i]*X[i,end] - (s[j])^2 for i in 1:N, j in 1:nieq) )
optimize!(model)
X1 = value.(X)
l  = @layout[a;b;c;d]
p11 = plot(LinRange(0,t[end],N), value.(v), title = "Bilevel case", label = ["v1" "v2" "v3" "v4"], xlabel = "t [h]", ylabel = "[mmol/(gdw*hr)]")
p12 = plot(LinRange(0,t[end],N), value.(z[:,:,1]), label = ["Glc" "Ac" "O2"], xlabel = "t [h]", ylabel = "[mmol]" )
p13 = plot(LinRange(0,t[end],N), value.(X[:,1]), label = "Biom.", xlabel = "t [h]", ylabel = "[g/l]")
p14 = plot(LinRange(0,t[end],N), value.(D), label = "Dilution rate", xlabel = "t [h]", ylabel = "[1/h]")
BiPlot = plot(p11,p12,p13,p14, layout = l)
return X1[:,1],z[:,:,1],v[:,:],D[:],s[:], BiPlot

elseif flag == 2
@NLobjective(model,Min,  + sum((μ[i]*L[neq+i+1]) for i in 1:nieq))
optimize!(model)
X2 = value.(X)
l  = @layout[a;b;c;d]
p21 = plot(LinRange(0,t[end],N), value.(v), title = "Regular case", label = ["v1" "v2" "v3" "v4"], xlabel = "t [h]", ylabel = "[mmol/(gdw*hr)]")
p22 = plot(LinRange(0,t[end],N), value.(z[:,:,1]), label = ["Glc" "Ac" "O2"], xlabel = "t [h]", ylabel = "[mmol]" )
p23 = plot(LinRange(0,t[end],N), value.(X[:,1]),  label = "Biom.", xlabel = "t [h]", ylabel = "[g/l]")
p24 = plot(LinRange(0,t[end],N), value.(D), label = "Dilution rate", xlabel = "t [h]", ylabel = "[1/h]")
BiPlot = plot(p21,p22,p23,p24,layout = l)
return X2[:,1],BiPlot

elseif flag == 3
@objective(model,Max, sum(D[i]*X[i,end] for i in 1:N))
optimize!(model)
X1 = value.(X)
p11 = plot(LinRange(0,t[end],N), value.(v), title = "Bilevel case", label = ["v1" "v2" "v3" "v4"], xlabel = "t [h]", ylabel = "[mmol/(gdw*hr)]")
p12 = plot(LinRange(0,t[end],N), value.(z[:,:,1]), label = ["Glc" "Ac" "O2"], xlabel = "t [h]", ylabel = "[mmol]" )
p13 = plot(LinRange(0,t[end],N), value.(X[:,1]), label = "Biom.", xlabel = "t [h]", ylabel = "[g/l]")
p14 = plot(LinRange(0,t[end],N), value.(D), label = "Dilution rate", xlabel = "t [h]", ylabel = "[1/h]")
@NLobjective(model,Min,  + sum((μ[i]*L[neq+i+1]) for i in 1:nieq))
optimize!(model)
X2 = value.(X)
p21 = plot(LinRange(0,t[end],N), value.(v), title = "Regular case", label = ["v1" "v2" "v3" "v4"], xlabel = "t [h]", ylabel = "[mmol/(gdw*hr)]")
p22 = plot(LinRange(0,t[end],N), value.(z[:,:,1]), label = ["Glc" "Ac" "O2"], xlabel = "t [h]", ylabel = "[mmol]" )
p23 = plot(LinRange(0,t[end],N), value.(X[:,1]),  label = "Biom.", xlabel = "t [h]", ylabel = "[g/l]")
p24 = plot(LinRange(0,t[end],N), value.(D), label = "Dilution rate", xlabel = "t [h]", ylabel = "[1/h]")
BiPlot = plot(p11,p21,p12,p22,p13,p23,p14,p24, layout = (4,2))
return X1[:,1],X2[:,1],BiPlot
end 
end