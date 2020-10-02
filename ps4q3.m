clear vars
close all
clc

% Parameters

N = 40;
mu = 0.5;
alpha = 0.4;
theta = 0.8;
delta = 0.3;
beta = 0.9;

% compute steady state 
K_ss = (beta*alpha*theta/(1 - beta*(1 - delta)))^(1/(1-alpha));
C_ss = theta*(K_ss^alpha) - delta*K_ss;

% Define grids for asset holdings
K = linspace(K_ss*0.5, K_ss*1.2, N)';

% Initialize two laws of motion
G1 = K_ss*ones(N,N);
G2 = K_ss*ones(N,N);
G1_new = zeros(N,N);
G2_new = zeros(N,N);

% compute income
K_grid = mu*repmat(K, 1, N) + (1-mu)*repmat(K', N, 1);
R = repmat(1 + (alpha*theta*(K_grid.^(alpha-1))) - delta, 1,1,N); 
W = repmat((1-alpha)*theta*(K_grid.^(alpha)),1,1,N);
Y = W + R.*repmat(reshape(K, 1, 1, N), N, N, 1);

% allocate space to make things faster
g = repmat(reshape(K, 1, 1, N), N, N, 1);
g_new = zeros(N,N,N);
C = ones(N,N,N)*C_ss;
C_new = zeros(N,N,N);
CP = ones(N,N,N);
RP = zeros(N,N,N);

% begin iteration
error = 1;
tolerance = 1e-4;
step = 0.2;
while error > tolerance
    
    pfi_error = 1;
    pfi_tolerance = 1e-4;
    step = 0.2;
        
    RP = 1 + alpha*theta*((mu*G1 + (1-mu)*G2).^(alpha-1)) - delta;
    
    % policy function iteration to find optimal g
    while pfi_error > pfi_tolerance
        % compute consumption and use Euler condition
        C = Y - g;
        CP = interp3(K, K, K, C, repmat(G1, 1,1,N), repmat(G2,1, 1, N), g, 'spline');
        C_new = (CP./repmat((beta*RP), 1, 1, N));
        % compute new policy g
        g_new = Y - C_new;
        % update g
        pfi_error = max(abs(g - g_new), [], 'all');
        g = (1-step)*g + step*g_new;
    end
        

    % Update the laws of motion
    for i = 1:N
        for j = 1:N
            G1_new(i,j) = g(i,j,i);
            G2_new(i,j) = g(i,j,j);
        end
    end
    error = max(max(abs(G1 - G1_new), [], 'all'), max(abs(G2 - G2_new), [], 'all'))
    G1 = (1-step)*G1 + step*G1_new;
    G2 = (1-step)*G2 + step*G2_new;
end

T = 200;
times = linspace(0, T, T);
K_A = zeros(1, T);
K_B = zeros(1, T);
K_C = zeros(1, T);

K_A(1) = 0.5*K_ss;
K_B(1) = 1.2*K_ss;
K_C(1) = mu*K_A(1) + (1-mu)*K_B(1);

for t = 2:T
    K_A(t) = interp3(K,K,K,g,K_A(t-1), K_B(t-1),K_A(t-1), 'spline');
    K_B(t) = interp3(K,K,K,g,K_A(t-1), K_B(t-1),K_B(t-1), 'spline');
    K_C(t) = mu*K_A(t) + (1-mu)*K_B(t);
end

f1 = figure;
plot(times, K_A)
hold on
plot(times, K_B)
hold on
plot(times, K_C)
saveas(f1, 'ps4q3.png')

