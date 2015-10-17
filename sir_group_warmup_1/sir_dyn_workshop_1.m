% George Zeng
% AM 115
% SIR Population Simulation

num_sim = 30;
num_indiv = 1000;
nc = 4;
pt = 0.5;
pr = 0.25;

% maximum reasnable number of weeks for calculating average
mat_size = 100;
avg_dynamics = zeros(mat_size,3);

figure();
hold on;
for i=1:num_sim
    [tt, results] = simulate1D(num_indiv,nc,pt,pr);
    plot(results);
    
    extra_zeros = size(avg_dynamics,1) - size(results,1);
    append_results = [results./num_sim ; zeros(extra_zeros,3)];
    avg_dynamics = avg_dynamics + append_results;
end

figure();
plot(avg_dynamics);
legend('Susceptible','Infected','Recovered');