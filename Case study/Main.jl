using JuMPDiff, Ipopt, Plots, Polynomials, Combinatorics
default(legend = :outertopright)
include("Colloc.jl")
include("MPC-Functions.jl")
include("DFBA-Functions.jl")
include("Bilevel-Functions.jl")

#Retrieve collocation matrix
adot  = colloc();
#Set up dictionary 
Setup = Dict(
    :t => 1:20,                     #Considered timespan
    :N => 100,                      #Number of finite elements
    :St => [0 -9.46 -9.48 -19.23;   #Glc     
            -39.34 0 1.24 12.12;    #Ac
            -35 -12.92 -12.73 0],   #O2
    :kla => 7.5,                    #[h^-1]
    :Km  => 0.015,                  #Michaelis-Menten constant
    :ncp => 3,                      #Default for colloc.jl, update here if changed
    :z0  => [10.8, 0.4, 0.21],      #Initial substrate concentrations
    :X0  => 0.001,                  #Initial biomass concentration
    :R   => 0.5,                    
    :Q   => 1
)


Setup;
X0  = Setup[:X0];
z0  = Setup[:z0];
N   = Setup[:N];
t   = Setup[:t];


#Comparison of different approaches for solving DFBA
#Note: Number of finite elements is redefined in DFBA_COMPARISON()
 DFBAplot = DFBA_COMPARISON(z0,X0)
 DFBAplot

#Bilevel optimization of DFBA in CSTR
#To do: Sensible penelization for slack variables
 X1,X2,BiPlot = Bilevel(z0,X0,3;N = 20, t = 10)
 BiPlot

#MPC
hor =  [3, 2]  #Prediction and control horizon
Din =  0;      #Dilution rate at t = 0
#Solve Bilevel optimization for setpoint 
Xsp,z,v,D,s,BiPlot = Bilevel(z0,X0,1,N = 40)
z,Xspsimp,X,Dm,v,MPCplot,MPCplot2 = MPCLoop(Din,Xsp,X0,z0,hor)
plot(LinRange(0,t[end],40),Xsp)
plot!(LinRange(0,t[end],N),Xspsimp)
plot!(LinRange(0,t[end],N),X)






