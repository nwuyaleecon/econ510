clc
clear
close all

% params
alpha = 0.66667;
beta = 0.98;
delta = 0.05;
A = 1;
theta = 2;
% steady state capital
k_s = (((1/beta)-(1-delta))/(alpha*A))^(1/(alpha-1));  
% initial capital
k_0 = 0.5*k_s;

% initialize all sequences
% length of sequence to approximate
T = 1000;
% capital sequence
k = linspace(k_0,k_s,T+2)';
% consumption sequence
c = zeros(T,1);
% capital investment return rate
r_k = zeros(T,1);
% wages
w = zeros(T+1,1);


% derived sequence
next_k = k;  

% algorithm execution parameters
max_error = 0.001;
max_iterations = 100000;
step = 0.0001;
error = 1.0;

% algorithm:
iteration_count = 0;
while ((iteration_count<max_iterations) && (error>max_error))
    
    % update consumption to be max feasible given capital (ensuring
    % positive consumption
    c(1:T+1)=max((A*k(1:T+1).^alpha)+(1-delta)*k(1:T+1)-k(2:T+2), 0);
    % update return rates according to euler condition
    r_k(1:T)=(((c(2:T+1).^theta)./c(1:T).^theta)/beta)-(1-delta);
    % update wages using firm maximization
    w(1:T+1)=k(1:T+1).^alpha*(1-alpha)*A;
    % use firm FOC
    next_k(2:T+1)=(r_k(1:T)./(alpha*A)).^(1/(alpha-1));
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
ylim(axes2,[0, 0.2]); 

figure3 = figure;
axes3 = axes('Parent',figure3);
plot(w)
xlim(axes3,[-10, T+10]);
ylim(axes3,[15, 35]); 