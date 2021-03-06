% George Zeng
% AM 115 
% Group Project 1
% 2D SIR

clc; clear all;

% Disease parameters
N = 200^2; % ideally a square number!
pI = 0.6;
pR = 0.4;
nsteps = 30;
% Population parameters
NI = 0.05*N; % number infected

% matrix of population individuals
WP = sqrt(N);
P = ones(WP, WP);

% Various infection initializations
%random points
I_ndx = randperm(N,NI);
P(I_ndx) = 2;
% square block
%P(:,40:60) = 2;

% Gaussian kernel - probability of infection decreases
% as distance to individual increases
sigma = 1;
NW = 10*sigma;
[Nx, Ny] = meshgrid(-NW:NW, -NW:NW);
G = 1./(sigma.*sqrt(2.*pi)).*exp(-(Nx.^2 + Ny.^2)./(2.*sigma^2));

S = double(P == 1);
I = double(P == 2);
R = double(P == 3);

% S_im = figure(1);
% I_im = figure(2);
% R_im = figure(3);

pop_im = figure(10);
pI_im = figure(11);

dyn_stats = [];

for i=1:nsteps
    
    % subpopulation indicator matrices
    S = double(P == 1);
    I = double(P == 2);
    R = double(P == 3);

    % show population dynamics with RGB color image
    P_im = zeros(WP,WP,3);
    P_im(:,:,1) = (S.*250 + I.*255 + R.*203);
    P_im(:,:,2) = (S.*250 + I.*0 + R.*203);
    P_im(:,:,3) = (S.*75  + I.*0 + R.*203);
    P_im = P_im./256;
    figure(pop_im);
    imshow(P_im);

    % stochastic component
    state = rand(WP,WP,1);
    pI_pop = pI.*conv2(I,G,'same');
    figure(pI_im); imshow(pI_pop); title('Probability of infection');
    
    S2 = S.*(state > pI_pop);
    I2 = S.*(state <= pI_pop) + I.*(state > pR);
    R2 = I.*(state <= pR) + R; 
    
    % individual dynamics figures
%     figure(S_im); imshow(S2); title('S');
%     figure(I_im); imshow(I2); title('I');
%     figure(R_im); imshow(R2); title('R');
    
    % new population matrix of individuals
    P = S2 + 2*I2 + 3*R2;
    
    % append current time period's stats
    i
    curr_stats = [sum(sum(S)), sum(sum(I)), sum(sum(R))];
    dyn_stats = [dyn_stats; curr_stats];
    
end

% show aggregate dynamics over time
dynamics_im = figure(20);
plot(dyn_stats);
legend('S','I','R');