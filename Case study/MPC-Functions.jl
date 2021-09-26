function BaseModel(Din,Xsp,z0,X0,hor,MPC;kwargs...)
    t   = get(kwargs, :t, Setup[:t]);
    N   = get(kwargs, :N, Setup[:N]);
    St  = get(kwargs, :St, Setup[:St]);
    kla = get(kwargs, :kla, Setup[:kla]);
    Km  = get(kwargs, :Km, Setup[:Km]);
    ncp = get(kwargs, :ncp, Setup[:ncp]);
    R   = get(kwargs, :R, Setup[:R]);
    Q   = get(kwargs, :Q, Setup[:Q]);
    ΔT  = t[end]/N;

    ## Set up the problem  
    np, nm = hor
    N = np;                          #Numer of finite elements equals horizon length
    ic   = 4;                        #Initial conditions
    neq  = 4*(N*ncp + (N-1)) + ic;   #Number of equality constraints
    nieq = 24*N;                     #Number of inequality constraints
    ncon = neq + nieq;               #Number of constraints
    zin = [10.8, 0.4, 0]

    
 
    
    model = Model(Ipopt.Optimizer);
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
    L[1] = sum(v[2,r] for i in 1:N,  r in 1:4) 
    #L[1]  = sum(X[i,r] for i in 1:N, r in 1:4) #Alternative formulation
    
    #Add intial conditions 
    for j in 1:3
        L[j+1] = (z0[j] - z[1,j,1]);   #Metabolite concentrations at t
        end 
    L[5] = (X0 - X[1,1]);              #Biomass concentration at t

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
    @variable(model,- 0.1 <= s[1:nieq] <= 0.01); #slack variables  
    @constraint(model, check[i in 22*N:23*N], s[i] >= -0.2);
    @constraint(model,check2[i in 23*N:24*N], s[i] >= -0.1)

    #KKT-constraints of the model
    @constraint(model,GRADL[i = 1:nvar-length(μ)-length(λ)], ∇L[i] == 0);
    @constraint(model,EQC[i in 1:neq], L[i+1] == 0);
    @NLconstraint(model,slack[i in 1:nieq], μ[i]*L[neq+1+i] + s[i]== 0 );
    @constraint(model, IEQC[i in neq+2:length(L)], L[i] >= 0);


    #MPC objective
    if MPC == true;
    @variable(model, dD[1:N], start = 0) 
    #MPC constraints 
    @constraint(model,dD[1] - (D[1] - Din) == 0);
    @constraint(model, [i in 2:N], dD[i] - (D[i] - D[i-1]) == 0);
    @NLobjective(model, Min, sum(ΔT*(0.5*(Q*(X[i,end]-Xsp)^2 + R*dD[i]^2)) + s[j]^2 for i=1:nm for j in 1:nieq));
    #@NLobjective(model,Min, sum(ΔT*(0.5*(Q*(X[i,end]-Xsp)^2 + R*dD[i]^2)) for i=1:nm) + 20*sum((μ[i]*L[neq+i+1]) for i in 1:nieq)) #Alternative formulation for weights in objective function
    optimize!(model)
    return  value(D[2]),value.(v[2,:])    #Optimal dilution and fluxes at next time step 
    elseif MPC == false;
    @constraint(model, [i in 1:N], Din == D[i]);
    @NLobjective(model,Min, + sum((μ[i]*L[neq+i+1]) for i in 1:nieq))
    optimize!(model)
    return value.(z[2,:,4]),value.(X[2,4]),value.(v[2,:])   #States at the next time step
    end 
end 

function MPCLoop(Din,Xsp,X0,z0,hor;kwargs...)
t   = get(kwargs, :t, Setup[:t]);
N   = get(kwargs, :N, Setup[:N]);
ΔT  = t[end]/N;
#State and control vectors
z = Array{Float64}(undef,(N,3));
X = Array{Float64}(undef,(N,1));
D = Array{Float64}(undef,(N,1));
v2 = Array{Float64}(undef,(N,4));
v1 = v2;
D[1]   = Din; 
z[1,:] = z0;
X[1]   = X0;
np,nc = hor;
#Convert the open loop optimal setpoint to a simpler setpoint with the correct Array length
top = findmax(Xsp[:,1]);
Xspsimp = Array{Float64}(undef,(N,1))
Xspsimp[1:oftype(1,round(top[2]*N/40))] .= X0; 
Xspsimp[oftype(1,round(top[2]*N/40)):end] .= 0.646; #TO DO: Set setpoint in more automated fashion - Problem: Max value not steady state value

for i in 1:N-1
    ## Optimal input using MPC 
    D[i+1],v1[i+1,:] = BaseModel(D[i],Xspsimp[i],z[i,:],X[i],[np,nc],true)
    ## Simulate plant with optimal input 
    z[i+1,:], X[i+1],v2[i+1,:] = BaseModel(D[i+1],Xspsimp[i],z[i,:],X[i],[np,0],false)
end 
p0 = plot(LinRange(0,t[end],N), value.(v1), label = ["v1" "v2" "v3" "v4"], xlabel = "t [h]", ylabel = "[mmol/(gdw*hr)]")
p1 = plot(LinRange(0,t[end],N), value.(v2), label = ["v1" "v2" "v3" "v4"], xlabel = "t [h]", ylabel = "[mmol/(gdw*hr)]")
p2 = plot(LinRange(0,t[end],N), value.(z[:,:,1]), label = ["Glc" "Ac" "O2"], xlabel = "t [h]", ylabel = "[mmol]" )
p3 = plot(LinRange(0,t[end],N), value.(X[:,1]), label = "Biom.", xlabel = "t [h]", ylabel = "[g/l]")
p4 = plot(LinRange(0,t[end],N), value.(D), label = "Dilution rate", xlabel = "t [h]", ylabel = "[1/h]")
p5 = plot(LinRange(0,t[end],N), value.(X[:,1]), label = "MPC", xlabel = "t [h]", ylabel = "[g/l]")
plot!(LinRange(0,t[end],N), value.(X[:,1]), label = "Open Loop")
MPCplot = plot(p0,p1,p2,p3,p4,p5,layout = (6,1))
return z,Xspsimp,X,D,v2,MPCplot,p5
end 
