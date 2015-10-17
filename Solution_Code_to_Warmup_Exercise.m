
% simulation of the stochastic SIR Model.

close all; clear all;

% set parameters for the SIR model
a = 2*log(2);
b = log(2);
Ntotal = 1000;

nsteps=20000;
numRuns=100;

%We solve the stochastic equation for 100 realization.

%the following variables will hold S(t) I(t) and R(t) for each realization in a matrix.
% the final size of these matrices will be nsteps*numRuns

SM=[];
IM=[];
RM=[];
figure(1)

for runs=1:numRuns
    runs
    %initialize vectors for S, I and R.
    S=zeros(nsteps,1);
    I=zeros(nsteps,1);
    R=zeros(nsteps,1);
    
    S(1)=Ntotal-10; % S populations size at the first time step is 990
    I(1)=10; % I populations size at the first time step is 10
    R(1)=0;
    
    t                   = 1; % initialize time
    dt                  = 1.e-3; % time step
    
    while (t<nsteps) % run loop and solve for each time step
        r=rand(1); % get a random number uniformly distributed between 0-1.
        
        % first case: One person got infected
        if r<((a/Ntotal)*I(t)*S(t)*dt)
            I(t+1)=I(t)+1;
            S(t+1)=S(t)-1;
            R(t+1)=R(t);
            
        % second case: One person recovered
        elseif r<( ((a/Ntotal)*I(t)*S(t)*dt) + b*I(t)*dt)
            I(t+1)= I(t)-1;
            R(t+1)= R(t)+1;
            S(t+1)=S(t);
            
        % third case: neither an infection nor recovery occurs 
        else
            I(t+1)=I(t);
            S(t+1)=S(t);
            R(t+1)=R(t);
        end
        
        t = t+1; % increment time
        
    end
    
    if(mod(runs,10)==0) %every 10 runs plot the S I R vectors as a function of time
        plot(dt*(1:1000:length(S)),S(1:1000:end),'ob')
        hold on
        plot(dt*(1:1000:length(S)),I(1:1000:end),'ok')
        plot(dt*(1:1000:length(S)),R(1:1000:end),'og')
        xlabel('time','fontsize',18)
        ylabel('S, I, R','fontsize',18)
        set(gca,'fontsize',18)
        h=legend('S','I','R');
        set(h,'fontsize',18);
    end
    
    % save the current vectors to the corresponding matrix
    SM=[SM,S(:)];
    IM=[IM,I(:)];
    RM=[RM,R(:)];
    
end

%Analytic solution the SIR model
[ t,P ] = solveSIR(a,b,Ntotal,Ntotal-10,20 );

% Plot the analytic solution on the same plot
plot(t(1,1:end),P(1,1:end),'b','linewidth',3);
plot(t(1,1:end),P(2,1:end),'k','linewidth',3);
plot(t(1,1:end),P(3,1:end),'g','linewidth',3);

% Compare the average of the realizations with the analytic solution
figure(2);
plot(dt*(1:nsteps),mean(SM,2),'ob');hold on
plot(dt*(1:nsteps),mean(IM,2),'ok');
plot(dt*(1:nsteps),mean(RM,2),'og');
plot(t(1,1:end),P(1,1:end),'b','linewidth',3);
plot(t(1,1:end),P(2,1:end),'k','linewidth',3);
plot(t(1,1:end),P(3,1:end),'g','linewidth',3);

xlabel('time','fontsize',18)
ylabel('S, I, R','fontsize',18)
set(gca,'fontsize',18)
h=legend('average S','average I','average R', 'theory S', 'theory I', 'theory r');
set(h,'fontsize',18);


%%
% Calculate the evolution of the PDFs as a function of time

figure(3)
linestylevector=['r';'b';'k';'c';'m'];
tvector=[1 5000 10000 15000 20000];

h1=subplot(3,1,1);
h2=subplot(3,1,2);
h3=subplot(3,1,3);
for tindex=1:5
    t=tvector(tindex);
    %calculate histogram for S
    [y,x]=hist(SM(t,:),1:10:1000);
    axes(h1)
    plot(x,y,linestylevector(tindex),'linewidth',3);
    set(gca,'fontsize',18)
    t=tvector(tindex);
    [y,x]=hist(IM(t,:),1:10:1000);
    axes(h2)
    plot(x,y,linestylevector(tindex),'linewidth',3);
    title('evolution of PDF for I','fontsize',18);
    ylabel('Population PDF','fontsize',18)
    set(gca,'fontsize',18)
     h=legend(['t=',num2str(t/1000)]);
     set(h,'fontsize',14);
    t=tvector(tindex);
    [y,x]=hist(RM(t,:),1:10:1000);
    axes(h3)
    plot(x,y,linestylevector(tindex),'linewidth',3);
    xlabel('Population size','fontsize',18)
    set(gca,'fontsize',18)
    title('evolution of PDF for R','fontsize',18);
    pause;
end
% set(gca,'fontsize',18)

%%

figure(4)
linestylevector=['r';'b';'k';'c';'m'];
tvector=[1 5000 10000 15000 20000];

subplot(3,1,1)   
for tindex=1:5
    t=tvector(tindex);
    %calculate histogram for S
    [y,x]=hist(SM(t,:),1:10:1000);
    plot(x,y,linestylevector(tindex),'linewidth',3);hold on;
end
title('Evolution of PDF for S','fontsize',18);
set(gca,'fontsize',18)
subplot(3,1,2)
for tindex=1:5
    t=tvector(tindex);
    [y,x]=hist(IM(t,:),1:10:1000);
    plot(x,y,linestylevector(tindex),'linewidth',3);hold on;
end
title('evolution of PDF for I','fontsize',18);
ylabel('Population PDF','fontsize',18)
set(gca,'fontsize',18)
h=legend('t=0','t=5','t=10', 't=15','t=20');
set(h,'fontsize',14);
subplot(3,1,3)
for tindex=1:5
    t=tvector(tindex);
    [y,x]=hist(RM(t,:),1:10:1000);
    plot(x,y,linestylevector(tindex),'linewidth',3);hold on;
end
xlabel('Population size','fontsize',18)
set(gca,'fontsize',18)
title('evolution of PDF for R','fontsize',18);

% set(gca,'fontsize',18)

