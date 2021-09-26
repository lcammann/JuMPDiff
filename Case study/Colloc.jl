#Create collocation matrix Adot 
#This section is just copy and pasted from Carols collocation code 
function colloc(ncp = 3, method = "Radau")
proots = Dict();
proots["Radau"]    = [[1.],
                  [0.333333; 1.],
                  [0.155051; 0.644049; 1.],
                  [0.088588; 0.409467; 0.787659; 1.],
                  [0.057104; 0.276843; 0.583590; 0.860240; 1.]];
                  
proots["Legendre"] = [[.5],
                  [0.211325; 0.788675],
                  [0.112702; 0.5; 0.887298],
                  [0.069432; 0.330009; 0.669991; 0.930568],
                  [0.046910; 0.230765; 0.5; 0.769235; 0.953090]];


tt = vcat([0.], proots[method][ncp]);
l = Array{Any,1}(undef, ncp+1);
dl = copy(l);
adot = zeros(ncp+1, ncp+1);
comb = collect(combinations(Array(1:ncp+1), ncp));
for i in 1:ncp+1
    k = setdiff(Array(1:ncp+1), comb[i])[1];
    l[k] = fromroots(tt[comb[i]])/prod((tt[k] - tt[comb[i]][j]) for j in 1:ncp);
    dl[k] = Polynomials.derivative(l[k]);
    for j in 1:ncp+1
        adot[k,j] = dl[k](tt[j]);
    end
end
return adot
end 
