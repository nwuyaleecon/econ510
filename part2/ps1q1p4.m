clc
clear
close all

% params
alpha = 0.66667;
beta = 0.98;
delta = 1;
A = 2;
A_0 = 1;
theta = 0.5;
% steady state capital
k_s = (((1/beta)-(1-delta))/(alpha*A))^(1/(alpha-1));  
% initial capital
k_0 = (((1/beta)-(1-delta))/(alpha*A_0))^(1/(alpha-1));

% initialize all sequences
% length of sequence to approximate
T = 50;
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

fprintf("Iterating for 2nd theta \n");
% change theta, initialize new sequences
theta = 1;
k2 = linspace(k_0,k_s,T+2)';
c2 = zeros(T,1);
r_k2 = zeros(T,1);
next_k = k2;  

% rerun the algorithm. reset error and iteration count
error = 1.0;
iteration_count = 0;
while ((iteration_count<max_iterations) && (error>max_error))
    
    % update consumption to be max feasible given capital (ensuring
    % positive consumption
    c2(1:T+1)=max((A*k2(1:T+1).^alpha)+(1-delta)*k2(1:T+1)-k2(2:T+2), 0);
    % update return rates according to euler condition
    r_k2(1:T)=(((c2(2:T+1).^theta)./c2(1:T).^theta)/beta)-(1-delta);
    % use firm FOC
    next_k(2:T+1)=(r_k2(1:T)./(alpha*A)).^(1/(alpha-1));
    error=max(abs(k2-next_k));
    
    % updating
    k2=(1-step)*k2+step*next_k;
    iteration_count = iteration_count+1; 
end

fprintf("Iterating for 3rd theta \n");

% change theta, initialize new sequences
theta = 2;
k3 = linspace(k_0,k_s,T+2)';
c3 = zeros(T,1);
r_k3 = zeros(T,1);
next_k = k3;  

% rerun the algorithm. reset error and iteration count
error = 1.0;
iteration_count = 0;
while ((iteration_count<max_iterations) && (error>max_error))
    
    % update consumption to be max feasible given capital (ensuring
    % positive consumption
    c3(1:T+1)=max((A*k3(1:T+1).^alpha)+(1-delta)*k3(1:T+1)-k3(2:T+2), 0);
    % update return rates according to euler condition
    r_k3(1:T)=(((c3(2:T+1).^theta)./c3(1:T).^theta)/beta)-(1-delta);
    % use firm FOC
    next_k(2:T+1)=(r_k3(1:T)./(alpha*A)).^(1/(alpha-1));
    error=max(abs(k3-next_k));
    
    % updating
    k3=(1-step)*k3+step*next_k;
    iteration_count = iteration_count+1; 
end

% plotting
figure1 = figure;
axes1 = axes('Parent',figure1);
plot(k);
hold on
plot(k2);
hold on
plot(k3);
xlim(axes1,[-10, T+10]);
ylim(axes1,[0, 1.5*k_s]); 
legend("\theta=0.5","\theta=1","\theta=2")
saveas(figure1, "ps1q1fig4.png")


s = (A*(k(1:T).^alpha) - c(1:T))./(A*(k(1:T).^alpha));
s2 = (A*(k2(1:T).^alpha) - c2(1:T))./(A*(k2(1:T).^alpha));
s3 = (A*(k3(1:T).^alpha) - c3(1:T))./(A*(k3(1:T).^alpha));

figure2 = figure;
axes2 = axes('Parent',figure2);
plot(s);
hold on
plot(s2);
hold on
plot(s3);
xlim(axes2,[-10, T+10]);
ylim(axes2,[0, 1]); 
legend("\theta=0.5","\theta=1","\theta=2")
saveas(figure2, "ps1q1fig5.png")
