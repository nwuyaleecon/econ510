clc 
clear
close all

% params
alpha = 0.3;
beta = 0.95;
theta = 1;

% compute steady states, just so that our grid occurs on interesting values
% of k, l
l_ss = (1 - alpha)/(2 - alpha - (alpha*beta));
k_ss = ((beta*alpha*theta)^(1/(1-alpha)))*l_ss;

% create the grids
k_grid_size = 500;
l_grid_size = 100;
k_grid = linspace(0.5*k_ss, 1.5*k_ss, k_grid_size);
l_grid = linspace(0.1, 0.9, l_grid_size);

% figure out what the optimal labor is at each k, k'
l_opt = zeros(k_grid_size, k_grid_size);
for i = 1:k_grid_size
    for j = 1:k_grid_size
        k = k_grid(i);
        k_prime = k_grid(j);
        possible_c = theta*(k^alpha)*(l_grid.^(1-alpha)) - k_prime;
        labor_objective = log(possible_c) + log(1-l_grid);
        [maxval, maxarg] = max(labor_objective, [], 'omitnan');
        l_opt(i,j) = maxarg;
    end
end

% figure out the utility grid at k, k'
U_grid = zeros(k_grid_size, k_grid_size);
for i = 1:k_grid_size
    for j = 1:k_grid_size
        k = k_grid(i);
        k_prime = k_grid(j);
        l = l_grid(l_opt(i, j));
        c = theta*(k^alpha)*(l^(1-alpha)) - k_prime;
        if c > 0
            U_grid(i,j) = log(c) + log(1-l);
        else
            U_grid(i,j) = -1e6;
        end
        
    end
end


% loop parameters
error = 1;
tolerance = 1e-8;

% value functions, correspondences
V_new = zeros(k_grid_size, 1);
V_old = zeros(k_grid_size, 1);
gk = zeros(k_grid_size, 1);

while error > tolerance
    [V_new, gk] = max(U_grid + beta*repmat(V_old', k_grid_size,1), [], 2);
    error = max(abs(V_old - V_new))
    V_old = V_new;
end

% check that this is close to analytical policy func
k_analytical = alpha*beta*theta*(k_grid.^alpha).*(l_ss.^(1-alpha));

f1 = figure;
plot(k_grid, k_grid(gk), k_grid, k_analytical);
xlabel('k')
ylabel('g(k)')

