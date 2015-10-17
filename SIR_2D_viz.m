% George Zeng
% AM 115 
% Group Project 1
% 2D dynamics of disease

clc; clear all;

% Disease parameters
pI = 0.6;
pR = 0.3;
nsteps = 30;
% Population parameters
WP = 100;
NI = 500;

% Initial infected
P = ones(WP, WP);
I_ndx = randperm(WP.^2,NI);
P(I_ndx) = 2;

% Gaussian kernel
sigma = 1;
NW = 10*sigma;
[Nx, Ny] = meshgrid(-NW:NW, -NW:NW);
G = 1./(sigma.*sqrt(2.*pi)).*exp(-(Nx.^2 + Ny.^2)./(2.*sigma^2));
%G(NW/2,NW/2) = 0;

S = double(P == 1);
I = double(P == 2);
R = double(P == 3);

S_im = figure(1);
I_im = figure(2);
R_im = figure(3);

pop_im = figure(10);

dyn_stats = [];

for i=1:nsteps
    
    % subpopulation indicator matrices
    S = double(P == 1);
    I = double(P == 2);
    R = double(P == 3);

    % show population dynamics;
    P_im = zeros(WP,WP,3);
    P_im(:,:,3) = S.*255;
    P_im(:,:,1) = I.*255;
    %P_im(:,:,2) = R.*255;
    figure(pop_im);
    imshow(P_im);

    % stochastic component
    state = rand(WP,WP,1);
    pI_pop = pI.*conv2(I,G,'same');
    %figure(pI_im); imshow(pI_pop); title('Probability of infection');
    
    S2 = S.*(state > pI_pop);
    I2 = S.*(state <= pI_pop) + I.*(state > pR);
    R2 = I.*(state <= pR) + R; 
    
    figure(S_im); imshow(S2); title('S');
    figure(I_im); imshow(I2); title('I');
    figure(R_im); imshow(R2); title('R');
    
    % new matrix
    P = S2 + 2*I2 + 3*R2;
    
    i
    curr_stats = [sum(sum(S)), sum(sum(I)), sum(sum(R))];
    dyn_stats = [dyn_stats; curr_stats];
    
end

dynamics_im = figure(20);
plot(dyn_stats);
legend('S','I','R');