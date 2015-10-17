
% simulation of the stochastic logistic equations.
% For APM 115, based on "modelling biological populations in space
% and time, Eric Renshaw.

close all; clear all

% set parameters for logistic equation: X_t=R*X*(1-X/K)
A                        = 2.*log(2);
B                        = log(2);
N                        = 1000;
dimP                     = N*10; % number of points at which pdf is defined.  Needs to be larger than carrying capacity:

%corresponding birth-death probability
%see lecture notes for the mapping
mu                       = B;
lambda                   = A/N;
 
% simulate the stochastic logistic equation, using section 3.4.
nsteps              = 4000;
I                   = zeros(nsteps, 1);
S                   = zeros(nsteps, 1);
R                   = zeros(nsteps, 1);

% variables to record probabilities of infection and recovery over time
pI                   = zeros(nsteps, 1);
pS                   = zeros(nsteps, 1);
pR                   = zeros(nsteps, 1);

% initial conditions
I(1)                = 1;
S(1)                = N - I(1);

t                   = 1;
dt                  = 1/100;
%Think about what happens if dt is too large. And why?
while (t<nsteps)
    rvI=rand(1);
    rvR=rand(1);
    pI(t) = lambda*I(t)*S(t)*dt; % infection probability
    pR(t) = mu*I(t)*dt; % recovery probability
    
    next_I = I(t);
    next_S = S(t);
    next_R = R(t);
  
  if (rvI < pI(t))
    next_I = I(t)+1; % corresponding to infection
    next_S = S(t)-1;
  end
  if (rvR < pR(t))
    next_I = I(t)-1; % corresponding to recovery
    next_R = R(t)+1;
  end
  
  I(t+1) = next_I;
  S(t+1) = next_S;
  R(t+1) = next_R;
  
  t = t+1;
end

figure();
hold on;
plot(I,'r-');
plot(S,'b-');
plot(R,'g-');
legend('Infected','Susceptible','Recovered');
title('SIR sub-population dynamics');
print -dpng 'gp_sir_stoch_pop.png';

figure();
hold on;
plot(pI,'r-');
plot(pR,'b-');
legend('p(Infection)','p(Recovery)');
title('SIR transition probabilities');
print -dpng 'gp_sir_stoch_probs.png';

%{
figure;
subplot(2, 1, 1)
plot(I(1:nsteps), '-r'); %an example of the time series
ylabel('I(t)');
xlabel('t');
title('stochastic simulation');

% calculate pdf from final, stationary, part of the time series:
pdf_simulation=zeros(dimP,1);
for t=floor(nsteps*3/4+eps):nsteps
  pdf_simulation(I(t))= pdf_simulation(I(t))+1;
end

subplot(2, 1, 2)
bh = bar(pdf_simulation./max(pdf_simulation), 'r');
set(bh, 'EdgeColor', 'k')
title('pdf from numerical simulation');
xlabel('N')
ylabel('number of occurences');
set(gca, 'XLim', [0, max(I)]);
pdf_simulation      = pdf_simulation/sum(pdf_simulation);
hold on;

% Calculate probability distribution from the analytic theory
NVec                     = 1:1:N;
pi_N                     = ((lambda/mu).^(NVec))./(factorial(NVec))./(exp(lambda/mu)-1);

pause
plot(NVec, pi_N./max(pi_N), '-b', 'LineWidth', 3);

%}