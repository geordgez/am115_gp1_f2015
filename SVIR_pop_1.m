% George Zeng
% AM 115
% Group Project
% Running SolveSVIR

N = 100000;
beta = 0.0000029; %2*log(2)/N;
gamma = 0.0026; %log(2);
e = 0.9;
p = 0.9;
mu = 0.01;
R0 = 1;
S0 = N-1;
tmax = 150;

[t,P]=solveSVIR(beta,gamma,e,p,mu,N,R0,S0,tmax);