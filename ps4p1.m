clc 
clear
close all

% params
alpha = 0.36;
beta_A = 0.95;
beta_B = 0.90;
theta = 1;
mu = 0.3;
delta = 0.05;
min_A = -2;
max_A = 4;

% using work from 1.2 to compute steady states. 
% Let the consumption ratio c_1 to c_2 to c, and denote the k1 to k2 ratio
% as r. let the capital-to-labor ratio be h, and let the l1 to l2 ratio be
% l. then we have
K_ss = (1/(alpha*theta*beta_A) - (1-delta)/(alpha*theta))^(1/(alpha-1));
A_s1 = (K_ss / mu) - min_A*((1-mu)/mu);
A_s2 = min_A;
W_ss = (1-alpha)*theta*K_ss^alpha;
R_ss = alpha*theta*K_ss^(alpha-1) + 1 - delta;

fprintf("Finished computing steady state \n")

% create the grid
A_size = 1000;
A_grid = linspace(min_A, max_A, A_size);

% populate interim utility
U_grid = zeros(A_size, A_size);

for i = 1:A_size
    val = log(max(R_ss*A_grid(i) + W_ss - A_grid, 1e-10));
    U_grid(i,:) = val ;
end

% loop parameters
error = 1;
tolerance = 1e-9;

% value functions
V_new = zeros(A_size, 1);
V_old = zeros(A_size, 1);
g_A = zeros(A_size, 1);
g_B = zeros(A_size, 1);


fprintf("Iterate to find g_A \n")

while error > tolerance
    for i = 1:A_size
        [V_new(i), g_A(i)] = max(U_grid(i,:)+ beta_A*V_old');
    end 
    error = max(max(abs(V_old - V_new)));
    V_old = V_new;
end

fprintf("Finished g_A, reiterate value function to find g_B \n")
% reset error
error = 1;
while error > tolerance
    for i = 1:A_size
        [V_new(i), g_B(i)] = max(U_grid(i,:)+ beta_B*V_old');
    end 
    error = max(max(abs(V_old - V_new)));
    V_old = V_new;
end

f1 = figure;
plot(A_grid,A_grid(g_A))
hold on
plot(A_grid,A_grid(g_B))
hold on
saveas(f1, 'ps4q1.png')
