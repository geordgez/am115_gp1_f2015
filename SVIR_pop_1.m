% George Zeng
% AM 115
% Group Project
% Running SolveSVIR

N = 100000;
beta = 2*log(2)/N;
gamma = log(2);
e = 0.8;
p = 0.91;
mu = 0.01;
R0 = 1;
S0 = N-1;
tmax = 100;

[t,P]=solveSVIR(beta,gamma,e,p,mu,N,R0,S0,tmax);