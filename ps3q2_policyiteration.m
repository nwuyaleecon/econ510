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
c_ss = theta*(k_ss^alpha)*(l_ss^(1-alpha)) - k_ss ;

% create the grids
k_grid_size = 20;
k_grid = linspace(0.8*k_ss, 1.2*k_ss, k_grid_size)';

% loop parameters
error = 1;
tolerance = 1e-8;
step = 0.3;
% suppress the output text on fsolve:
options = optimoptions('fsolve','Display','off');

% policy functions
c_k = ones(k_grid_size, 1)*c_ss;
c_interpolated = zeros(k_grid_size, 1);
k_prime = zeros(k_grid_size, 1);
l_k = zeros(k_grid_size, 1);
l_prime = zeros(k_grid_size, 1);


while error > tolerance
    % for each capital level in grid k, find optimal labor l
    
    for i = 1:k_grid_size
        k = k_grid(i);
        c = c_k(i);
        l_k(i) = fsolve(@(l)optimal_labor_function(l,k,c), l_ss, options);
        k_prime(i) = theta*(k^alpha)*(l_k(i)^(1-alpha)) - c;
        interpolated_value = interp1(k_grid, c_k, k_prime(i), 'linear', 'extrap');
        c_interpolated(i) = interpolated_value;
        k = k_prime(i);
        c = c_interpolated(i);
        l_prime(i) = fsolve(@(l)optimal_labor_function(l,k,c), l_ss, options);
    end
    
    c_new = c_interpolated ./ (beta * alpha * theta * (k_prime.^(alpha-1)).*(l_prime.^(1-alpha)));
    error = max(abs(c_new - c_k))
    c_k = step*c_new + (1-step)*c_k;
end



% check that this is close to analytical policy func
k_analytical = alpha*beta*theta*(k_grid.^alpha).*(l_ss.^(1-alpha));

f1 = figure;
plot(k_grid, k_prime, k_grid, k_analytical);
xlabel('k')
ylabel('g(k)')

f2 = figure;
plot(k_grid, (k_prime - k_analytical))
xlabel('k')
ylabel('error')

function opt_labor = optimal_labor_function(l, k, c)
    alpha = 0.3;
    theta = 1;
    opt_labor = (1-alpha)*theta*(k^alpha)*(l^(-1*alpha))*(1-l) - c;
end

