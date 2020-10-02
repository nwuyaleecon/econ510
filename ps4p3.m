

% params
alpha = 0.4;
beta = 0.9;
theta = 0.8;
mu = 0.5;
delta = 0.3;
% using work from 1.2 to compute steady states. 
% Let the consumption ratio c_1 to c_2 to c, and denote the k1 to k2 ratio
% as r. let the capital-to-labor ratio be h, and let the l1 to l2 ratio be
% l. then we have
K_ss = ((alpha*theta*beta)/(1 - beta*(1 - delta)))^(1/(1-alpha));
C_ss = theta*K_ss^alpha-delta*K_ss;
min_K = K_ss*0.5;
max_K = K_ss*1.2;


fprintf("Finished computing steady state \n")

% create the grid
K_size = 50;
K_grid = linspace(min_K, max_K, K_size);
R_grid = alpha*theta*(K_grid.^(alpha - 1)) + (1 - delta);
W_grid = (1-alpha)*theta*(K_grid.^(alpha));

% at each level of aggregate capital, this is the income
Y = repmat(R_grid, K_size, 1).*repmat(K_grid', 1, K_size) + repmat(W_grid, K_size, 1); 
C = zeros(K_size, K_size);

% loop parameters
error = 1;
tolerance = 1e-5;
pfi_tolerance = 1e-5;
step_G = 0.2;

% laws of motion
G = ones(K_size, 1)*K_ss;
G_new = ones(K_size, 1)*K_ss;

% consumption policy functions
c_old = ones(K_size, K_size)*C_ss;
c_new = ones(K_size, K_size)*C_ss;

% capital policy function
g = ones(K_size, K_size);

fprintf("Iterate to find g \n")
while error > tolerance
    
    % iterate over the value function first to find the fixed point of the
    % contraction
    pfi_error = 1;
    step = 0.2;
    
    
    
    while pfi_error > pfi_tolerance
        for i=1:K_size
                k = K_grid(i);
                k_prime=(theta*alpha*G.^(alpha-1)+(1-delta))*k + (1-alpha)*(G.^alpha)- c_old(i,:)';        
                c_prime=interp2(K_grid,K_grid,c_old,k_prime,G,'spline');        
                c_new(i,:)=c_prime/beta./(theta*alpha*k_prime.^(alpha-1)+1-delta);
            
        end
    %    C = max(Y - g_old, 1e-6);
    %    g_interp = interp2(K_grid, K_grid, g_old, g_old, repmat(G,1,K_size), 'spline');
    %    
    %    R_prime = 1 + alpha*theta*(G.^(alpha-1)) - delta;
    %    W_prime = (1-alpha)*(G.^alpha);
    %    
    %    Y_prime = repmat(W_prime, 1, K_size) + repmat(R_prime, 1, K_size) .* g_old;
    %    C_prime = max(Y_prime - g_interp, 1e-6);
    %    
    %    C_new = C_prime./(beta*R_prime);
    %    g_new = Y - C_new;
        pfi_error = max(abs(c_old - c_new), [], 'all');
        c_old = (1-step)*c_old + step*c_new;
    end
    
    for i = 1:K_size
        for j = 1:K_size
            K = K_grid(j);
            k = K_grid(i);
            g(i,j) = k*(alpha*theta*K^(alpha-1) + 1 - delta)+ (1-alpha)*theta*K^alpha - c_old(i,j);
        end
    end
    
    % now adjust G
    for K = 1:K_size
        G_new(K) = g(K, K);
    end
    
    error = max(max(abs(G_new - G)))
    G = step_G * G_new + (1-step_G)*G;
    
end

g_policy = ones(K_size, K_size);

for i = 1:K_size
    for j = 1:K_size
        k = K_grid(i);
        K = K_grid(j);
        g_policy(i,j) = k*(theta*alpha*K^(alpha - 1) + 1 - delta) + (1-alpha)*theta*K^alpha- c_old(i,j);
    end
end

% we have the aggregate law of motion on G, and g_old is the optimal policy
% function consistent with this law of motion

f1 = figure;
plot(K_grid,G)
%saveas(f1, 'ps4q1.png')

T = 200;
times = linspace(1, T, T);
K_A = zeros(T, 1);
K_B = zeros(T, 1);
K_C = zeros(T, 1);

K_A(1) = 0.5*K_ss;
K_B(1) = 1.2*K_ss;
K_C(1) = mu*K_A(1) + (1-mu)*K_B(1);

for t = 2:T
    prev_K = mu*K_A(t-1)+(1-mu)*K_B(t-1);
    K_A(t) = interp2(K_grid, K_grid, g_policy, prev_K, K_A(t-1), 'spline');
    K_B(t) = interp2(K_grid, K_grid, g_policy, prev_K, K_B(t-1), 'spline');
    K_C(t) = interp1(K_grid, G, K_C(t-1));
end

f1 = figure;
plot(times,K_A)
hold on
plot(times,K_B)
hold on
plot(times, K_C)
saveas(f1, 'ps4q3.png')

f2 = figure;
surf(K_grid, K_grid, g_policy)


