clc
clear
close all

r = 0.02;
beta = 0.96;
min_A = -0.1;
epsilon = 0.1;
p = 0.95;

% initialize grid of capital
K_min = min_A;
K_max = 1;
N = 200;
K = linspace(K_min, K_max, N);

% precompute U at each pair k, k'
fprintf("precomputing utilities \n")
U_H = zeros(N,N);
U_L = zeros(N,N);
for i = 1:N
    for j = 1:N
        k = K(i);
        kp = K(j);
        U_H(i,j) = log(max(1 + epsilon + (1+r)*k - kp, 1e-10));
        U_L(i,j) = log(max(1 - epsilon + (1+r)*k - kp, 1e-10));
    end
end

fprintf("value function iteration \n")
V_H = zeros(N,1);
V_L = zeros(N,1);
V_H_new = zeros(N,1);
V_L_new = zeros(N,1);
error = 1;
tolerance = 1e-6;

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


fprintf("value function iteration, varying min holding  \n")
min_A = -0.2;
K_min = min_A;
K_max = 1;
N = 200;
K2 = linspace(K_min, K_max, N);

U_H = zeros(N,N);
U_L = zeros(N,N);
for i = 1:N
    for j = 1:N
        k = K2(i);
        kp = K2(j);
        U_H(i,j) = log(max(1 + epsilon + (1+r)*k - kp, 1e-10));
        U_L(i,j) = log(max(1 - epsilon + (1+r)*k - kp, 1e-10));
    end
end

V_H = zeros(N,1);
V_L = zeros(N,1);
V_H_new = zeros(N,1);
V_L_new = zeros(N,1);
error = 1;
tolerance = 1e-6;

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


f1 = figure;
plot(K, K(g1_H))
hold on
plot(K, K(g1_L))
hold on
plot(K2, K2(g2_H))
hold on
plot(K2, K2(g2_L))
hold on
h = refline(1,0);
h.Color = 'black';
legend({'g_H, bmin=-0.1', 'g_L, bmin=-0.1','g_H, bmin=-0.2', 'g_L, bmin=-0.2', '45 deg line' }, 'Location', 'southeast')
saveas(f1, 'ps5q3p3_3.png')