function [ t,P ] = solveSIR(a,b,N,S0,tmax )
%N is total population
%S0 is the initial suseptible population
%
%in the logistic model we had only a single dependent
%variable. here we have 3 dependent variables and we will
%refer to them as a vector P=[S,I,R]
%the individual variables are P(1)=1, P(2)=I and P(3)=R

sol=ode45(@SIRmodel,[0 tmax],[S0,1,N-S0-1]);
%the 3 inital values are given as a vector to ode45
%we are assuming that at t=0 1 individual is infected and
%all the others are susceptible
%SIRmodel is defined below

t=linspace(0,tmax); %100 equally spaced points
P=deval(sol,t);

%now we can plot the solution
plot(t,P)
%remember P(t) is a vector of 3 values so this plot will have
%all 3 cureves on it
%we can select individual variables from the t and P returned
%by solveSIR in the command window

%ALWAYS label your axes
xlabel('time')
ylabel('numbers of individuals')
%a title wouldn't hurt either
title('S,I and R vs time, t')

%and a legend allows us to tell the curves apart
legend('S','I','R')
%produces a legend associating the given strings with the
%3 curves in their order in the P vector

    function dP=SIRmodel(t,P)
        %remember P is a vector [S,I,R] and we must calculate
        %the derivatives as a column vector [dS/dt;dI/dt;dR/dt]
        dP=zeros(3,1);
        %a column vector of zeros (matrix of 3 rows, 1 column)
        %into which we place our derivatives
        dP(1)=-a*P(1)*P(2)/N;
        dP(2)=a*P(1)*P(2)/N - b*P(2);
        dP(3)=b*P(2);
    end   %ends our derivatives function
end   %ends our main solver function

%with this solver function saved in our current directory as a .m file
%we could then do a solution for particular values of the
%parameters in the command window as
%  [t,P]=solveSIR(2*0.693,0.693,1000,999,30);
