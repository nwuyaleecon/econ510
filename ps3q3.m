clc 
clear
close all

% params
alpha = 0.36;
beta = 0.95;
theta = 1;
sigma = 2;
mu = 0.5;
delta_1 = 0.05;
delta_2 = 0.08;


% using work from 3.3 compute steady states. 
% Let the consumption ratio c_1 to c_2 to c, and denote the k1 to k2 ratio
% as r. let the capital-to-labor ratio be h, and let the l1 to l2 ratio be
% l. then we have
h_1 = (beta*theta*alpha/(1 - beta*(1 - delta_1)))^(1/(1-alpha));
h_2 = (beta*theta*mu/(1 - beta*(1 - delta_2)))^(1/(1-mu));
c = (((1-alpha)*h_1^alpha)/((1-mu)*h_2^mu))^(1/sigma);
r = c*(h_2^(mu - 1) - delta_2)/(h_1^(alpha - 1) - delta_1);
l = r*h_2/h_1;
l_1 = 1/(1 + 1/l);
l_2 = 1 - l_1;
k_s1 = h_1 * l_1
k_s2 = h_2 * l_2



fprintf("Finished computing steady state \n")

% create the grids
k_grid_size = 75;
l_grid_size = 75;
k_min = min(k_s1, k_s2);
k_max = max(k_s1, k_s2);
k_grid = linspace(0.4*k_min, k_max, k_grid_size);
l_grid = linspace(0.3, 0.7, l_grid_size);


fprintf("Populating optimal labor choices, and subsequent utilities \n")

% figure out what the optimal labor is at each k_1, k_2, k_1', k_2'
l_opt = zeros(k_grid_size, k_grid_size, k_grid_size, k_grid_size);
U_grid = zeros(k_grid_size, k_grid_size, k_grid_size, k_grid_size);
y_1 = zeros(l_grid_size, 1);
y_2 = zeros(l_grid_size, 1);
c_1 = zeros(l_grid_size, k_grid_size, k_grid_size);
c_2 = zeros(l_grid_size, k_grid_size, k_grid_size);
c_min = 1d-8*ones(l_grid_size, k_grid_size, k_grid_size);

for i1 = 1:k_grid_size
    for i2 = 1:k_grid_size
                
        k1 = k_grid(i1);
        k2 = k_grid(i2);
        y_1 = theta*(k1^alpha)*(l_grid.^(1-alpha)) + (1-delta_1)*k1;
        y_2 = theta*(k2^mu)*((1-l_grid).^(1-mu))+ (1-delta_2)*k2;
        
        c_1 = repmat(repmat(y_1', 1, k_grid_size) - repmat(k_grid,l_grid_size,1), 1, 1, k_grid_size);
        c_1 = max(c_1, c_min);
        
        c_2 = repmat(reshape(repmat(y_2', 1, k_grid_size) - repmat(k_grid,l_grid_size,1), l_grid_size, 1, k_grid_size), 1, k_grid_size, 1);
        c_2 = max(c_2, c_min);
        
        [U_grid(i1,i2,:,:), l_opt(i1,i2,:,:)] = max(((c_1.^(1-sigma))./(1-sigma)) + ((c_2.^(1-sigma))./(1-sigma)), [], 1);
    end
end

% loop parameters
error = 1;
tolerance = 1e-6;

% value functions, correspondences
V_new = zeros(k_grid_size, k_grid_size);
V_old = zeros(k_grid_size, k_grid_size);

fprintf("Starting iteration \n")

while error > tolerance
    for i = 1:k_grid_size
        for j = 1:k_grid_size
            V_new(i,j) = max(reshape(U_grid(i,j,:,:), [k_grid_size,k_grid_size])+ beta*V_old, [], 'all');
        end
    end
    error = max(max(abs(V_old - V_new)))
    V_old = V_new;
end

% find policy functions
g1 = zeros(k_grid_size, k_grid_size);
g2 = zeros(k_grid_size, k_grid_size);
gL = zeros(k_grid_size, k_grid_size);

for i = 1:k_grid_size
    for j = 1:k_grid_size
        
        V = max(reshape(U_grid(i,j,:,:), [k_grid_size,k_grid_size]) + beta*V_old, [], 'all');
        [g1(i,j), g2(i,j)] = find(reshape(U_grid(i,j,:,:), [k_grid_size,k_grid_size]) + beta*V_old == V);
        gL(i,j) = l_grid(l_opt(i,j,g1(i,j), g2(i,j)));
        
    end
end

f1 = figure;
surf(k_grid, k_grid, V_old)

f2 = figure;
surf(k_grid, k_grid, k_grid(g1))

f3 = figure;
surf(k_grid, k_grid, k_grid(g2))

f4 = figure;
surf(k_grid, k_grid, gL)

% now plot sequences of capital, consumption, labor
T = 200;
K_1 = zeros(T, 1);
idx_1 = zeros(T,1);
K_2 = zeros(T, 1);
idx_2 = zeros(T, 1);
C_1 = zeros(T, 1);
C_2 = zeros(T, 2);
times = linspace(1, T, T);
L = zeros(T, 1);

% populate initial values
[~, idx_1(1)] = min(abs(k_grid - 0.5*k_s1));
[~, idx_2(1)] = min(abs(k_grid - 0.5*k_s2));
K_1(1) = k_grid(idx_1(1));
K_2(1) = k_grid(idx_2(1));
L(1) = gL(idx_1(1), idx_2(1));

for t = 2:200
    idx_1(t) = g1(idx_1(t-1), idx_2(t-1));
    idx_2(t) = g2(idx_1(t-1), idx_2(t-1));
    K_1(t) = k_grid(idx_1(t));
    K_2(t) = k_grid(idx_2(t));
    L(t) = gL(idx_1(t), idx_2(t));
    C_1(t-1) = (K_1(t-1)^alpha)*(L(t-1)^(1-alpha)) + (1-delta_1)*K_1(t-1) - K_1(t);
    C_2(t-1) = (K_2(t-1)^mu)*(L(t-1)^(1-mu)) + (1-delta_2)*K_2(t-1) - K_2(t);
end

f5 = figure;
plot(times,K_1)
hold on
plot(times,K_2)
saveas(f5, 'ps3q3_k.png')

f6 = figure;
plot(times,L)
hold on
plot(times,1-L)
saveas(f6, 'ps3q3_l.png')

f7 = figure;
plot(times,C_1)
hold on
plot(times,C_2)
saveas(f7, 'ps3q3_c.png')
