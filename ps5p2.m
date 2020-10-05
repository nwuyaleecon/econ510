clc
clear
close all

% set constants
alpha = 0.3;
delta =  0.1;
beta = 0.96;
epsilon = 0.05;
p = 0.5;

% set capital grid
k_min = 0.9*((alpha*(1-eps))/(beta^(-1) - 1 + delta))^(1/(1-alpha));
k_max = 1.1*((alpha*(1+eps))/(beta^(-1) - 1 + delta))^(1/(1-alpha));
N = 500;
K = linspace(k_min, k_max, N);

% pre compute consumption and log-utility at each pair of k, k'
C_H = zeros(N, N);
C_L = zeros(N, N);
U_H = zeros(N, N);
U_L = zeros(N, N);

fprintf("precomputing utility/consumption matrix \n");

for i = 1:N 
    for j = 1:N
        k = K(i);
        kp = K(j);
        C_H(i,j) = (1+epsilon)*(k^alpha) + (1-delta)*k - kp;
        C_L(i,j) = (1-epsilon)*(k^alpha) + (1-delta)*k - kp;
        U_H(i,j) = log(max( C_H(i,j), 1e-10 ));
        U_L(i,j) = log(max( C_L(i,j), 1e-10 ));
    end
end


fprintf("starting value function iteration, p=0.5 \n");
% define two value functions, we iterate over the pair
V_H = zeros(N, 1);
V_L = zeros(N, 1);
V_H_new = zeros(N, 1);
V_L_new = zeros(N, 1);
error = 1;
tolerance = 1e-5;

while error > tolerance
    V_H_new = max(U_H + beta*p*repmat(V_H',N,1) + beta*(1-p)*repmat(V_L',N,1), [], 2);
    V_L_new = max(U_L + beta*p*repmat(V_L',N,1) + beta*(1-p)*repmat(V_H',N,1), [], 2);
    
    error = max(max(abs(V_L_new - V_L), [],'all'), max(abs(V_H_new - V_H), [], 'all'));
    V_L = V_L_new;
    V_H = V_H_new;
end

% extract policy functions
[~, g1_H] = max(U_H + beta*p*repmat(V_H',N,1) + beta*(1-p)*repmat(V_L',N,1), [], 2);
[~, g1_L] = max(U_L + beta*p*repmat(V_L',N,1) + beta*(1-p)*repmat(V_H',N,1), [], 2);

fprintf("starting value function iteration, p=0.95 \n");
% define two value functions, we iterate over the pair
p = 0.95;
V_H = zeros(N, 1);
V_L = zeros(N, 1);
V_H_new = zeros(N, 1);
V_L_new = zeros(N, 1);
error = 1;
tolerance = 1e-5;

while error > tolerance
    V_H_new = max(U_H + beta*p*repmat(V_H',N,1) + beta*(1-p)*repmat(V_L',N,1), [], 2);
    V_L_new = max(U_L + beta*p*repmat(V_L',N,1) + beta*(1-p)*repmat(V_H',N,1), [], 2);
    
    error = max(max(abs(V_L_new - V_L), [],'all'), max(abs(V_H_new - V_H), [], 'all'));
    V_L = V_L_new;
    V_H = V_H_new;
end

% extract policy functions
[~, g2_H] = max(U_H + beta*p*repmat(V_H',N,1) + beta*(1-p)*repmat(V_L',N,1), [], 2);
[~, g2_L] = max(U_L + beta*p*repmat(V_L',N,1) + beta*(1-p)*repmat(V_H',N,1), [], 2);

% Plot
f1 = figure;
plot(K, K)
hold on
plot(K, K(g1_H))
hold on
plot(K, K(g1_L))
hold on
plot(K, K(g2_H))
hold on
plot(K, K(g2_L))
legend({'45 degree line', 'z_H, p=0.5','z_L, p=0.5', 'z_H, p=0.95', 'z_L, p=0.95'},'Location','southeast')
saveas(f1, 'ps5q2.png')

