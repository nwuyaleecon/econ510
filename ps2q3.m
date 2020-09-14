clc
clear
close all

% params
alpha = 0.36;
beta = 0.9494;
delta = 0.0667;
theta = 0.6733;
sigma = 2;
% steady state capital
k_s = (((1/beta)-(1-delta))/(alpha*theta))^(1/(alpha-1));  
% initial capital
k_0 = 1.5*k_s;

% initialize all sequences
% length of sequence to approximate
T = 100;
% capital sequence
k = linspace(k_0,k_s,T+2)';
% consumption sequence
c = zeros(T,1);
% capital investment return rate
r_k = zeros(T,1);
% bond investment return rate
r_b = zeros(T,1);
% wages
w = zeros(T+1,1);


% derived sequence
next_k = k;  

% algorithm execution parameters
max_error = 0.001;
max_iterations = 5000;
step = 0.001;
error = 1.0;

% algorithm:
iteration_count = 0;
while ((iteration_count<max_iterations) && (error>max_error))
    
    % update consumption to be max feasible given capital
    c(1:T+1)=max((theta*k(1:T+1).^alpha)+(1-delta)*k(1:T+1)-k(2:T+2), 0);
    % update return rates according to euler condition
    r_k(1:T)=(((c(2:T+1).^sigma)./c(1:T).^sigma)/beta)-(1-delta);
    % update bond interest rate using no-arb
    r_b(1:T)=r_k(1:T) - delta;
    % update wages using firm maximization
    w(1:T+1)=k(1:T+1).^alpha*(1-alpha)*theta;
    % use firm FOC
    next_k(2:T+1)=(r_k(1:T)./(alpha*theta)).^(1/(alpha-1));
    error=max(abs(k-next_k));
    
    % updating
    k=(1-step)*k+step*next_k;
    iteration_count = iteration_count+1; 
end

% plotting
figure1 = figure;
axes1 = axes('Parent',figure1);
plot(k)
xlim(axes1,[-10, T+10]);
ylim(axes1,[0, 1.5*k_s]); 

figure2 = figure;
axes2 = axes('Parent',figure2);
plot(r_k)
xlim(axes2,[-10, T+10]);
ylim(axes2,[0, 0.5]); 

figure3 = figure;
axes3 = axes('Parent',figure3);
plot(r_b)
xlim(axes3,[-10, T+10]);
ylim(axes3,[0, 0.5]); 


figure4 = figure;
axes4 = axes('Parent',figure4);
plot(w)
xlim(axes4,[-10, T+10]);
ylim(axes4,[0, 1]); 