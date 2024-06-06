# ============================================================ #
#  JULIA CODE WRITTEN FOR CH7 OF KEIZAI SAMINAR                #
#  WRITTEN BY Sagiri Kitao (Translated in Julia by Taiki Ono)  #
#  09/13/2022                                                  #
# ============================================================ #

using NLsolve
using Plots
using LaTeXStrings
include("policyA.jl")

# create constructer that contains parameters
struct Model{TI<:Integer, TF<:AbstractFloat}
    
    alpha::TF      
    del::TF          
    beta::TF           
    Nj::TI        
    Njw::TI                
    rho::TF              
    meaJ::Array{TF,1}
    theta::Array{TF,1}
    L::TF
    ss0::Array{TF,1}
    tau::TF

end


ind_E = 4;
    # =1 BASELINE
    # =2 rho 50% down
    # =3 raise NRA TO 70
    # =4 live up to 85


# parameters
alpha = 0.40;
del = 0.08;
beta = 0.98;
Nj = 61;  # live 60 years max
Njw = 45; # working years (retires at Njw+1) : enter at 21 and retire at 65
rho = 0.5 # SS replacement rate (0.0 or 0.5)

meaJ = (1/Nj) .* ones(Nj); # age distribution (no population growth)

if ind_E == 2
    rho = 0.25;
end

if ind_E == 3
    Njw = 50;
end

if ind_E == 4
    Nj = 66;
    meaJ = (1/Nj) .* ones(Nj); # in experiment to raise longevity
end

theta = zeros(Nj);
theta[1:Njw] .= 1.0;

L = sum(meaJ[1:Njw] .* theta[1:Njw]);

ss0 = zeros(Nj);
ss0[Njw+1:Nj] .= rho;

tau = rho*(Nj-Njw)/Njw;


#  create constructer that contains parameters
m = Model(alpha,del,beta,Nj,Njw,rho,meaJ,theta,L,ss0,tau);

# initial guess for asset decisions
a_g = zeros(Nj-1); # a[2] to a[Nj] (NOTE: a[1]=a[Nj+1]=0) 
a_g[1] = 0.01;
for j in 1:m.Nj-2
    a_g[j+1] = a_g[j] + 0.01;
end

# find solutions: aX[2:Nj]
aX = nlsolve(x->policyA(m,x),a_g).zero;

a = zeros(m.Nj+1);
a[2:m.Nj] .= copy(aX);

K = sum(a[1:end-1] .* m.meaJ);

r = m.alpha * (K/L)^(m.alpha-1) - m.del; # interest rate
w = (1-m.alpha) * (K/L)^m.alpha; # wage

ss = ss0 .* w;

c = zeros(m.Nj);
for j in 1:m.Nj
    c[j] = w*m.theta[j]*(1-m.tau) + (1+r)*a[j] + ss[j] - a[j+1];
end

ageA = collect(1:Nj+1);
ageA .+= 19;

ageC = collect(1:Nj);
ageC .+= 19;

minJ = 20;
maxJ = 19 + m.Nj;

norm = c[1];

plot(ageA,a/norm,lw=3,legend=false)
title!(L"Asset: $a_{j}$")
xlims!(minJ,maxJ)
savefig("fig_a.pdf")

plot(ageC,c[1:m.Nj]./norm,lw=3,legend=false)
title!(L"Comsumption: $c_{j}$")
xlims!(minJ,maxJ)
savefig("fig_c.pdf")



