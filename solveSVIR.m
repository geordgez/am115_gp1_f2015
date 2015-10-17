function [ t,P ] = solveSVIR(beta,gamma,e,p,mu,N,R0,S0,tmax )

sol=ode45(@SIRmodel,[0 tmax],[S0,R0,1,N-S0-1]);

t=linspace(0,tmax); %100 equally spaced points
P=deval(sol,t);

%now we can plot the solution
plot(t,P)

%ALWAYS label your axes
xlabel('time')
ylabel('numbers of individuals')
%a title wouldn't hurt either
title('S, V, I and R vs time, t')

%and a legend allows us to tell the curves apart
legend('S','V','I','R')
%produces a legend associating the given strings with the
%3 curves in their order in the P vector

    function dP=SIRmodel(t,P)
        dP=zeros(4,1);
        
        dP(1)= (1-e*p)*mu*N - beta*P(1)*P(3) - mu*P(1); %S
        dP(2)= e*p*mu*N - mu*dP(2); % added V
        dP(3)= beta*P(1)*P(3) - gamma*P(3) - mu*P(3); %I
        dP(4)= gamma*P(3) - mu*P(4); %R
        
    end   %ends our derivatives function
end   %ends our main solver function

%with this solver function saved in our current directory as a .m file
%we could then do a solution for particular values of the
%parameters in the command window as
%  [t,P]=solveSIR(2*0.693,0.693,1000,999,30);
