% George Zeng
% AM 115
% Group Project
% Running SolveSVIR

N = 100000;
beta = 2*log(2)/N;
gamma = log(2);
e = 0.9;
p = 0.91;
mu = 0.02;
S0 = 99999;
tmax = 80;

[t,P]=solveSVIR(beta,gamma,e,p,mu,N,S0,tmax);