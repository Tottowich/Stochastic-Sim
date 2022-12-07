%% Task2: Back To Origin. (Brownian Motion)

%{      
 y=1|-----------|
  | |           |
  2 |     o     |
  | |           |
y=-1|-----------|
   x=-1 --2--  x=1
o: is the origin of the system (0,0)
This is an illustration of the box
%}

% A partical moves inside the box as: (x(t),y(t))=(W_1(t),W_2(t)), t>=0
% In this case we have standard BM=>sigma_1=sigma_2=1

%% Task Visualization A)
dt = 0.0001; % Specified in the assignment.
sigma = 1;
W = @(p) p+sqrt(dt)*sigma*randn(size(p)); % Function to update position.
boundaries = [-1,-1;1,1];
N = 1;
history = [];
for i=1:N
    xy = zeros(1,2); % Each simulation starts at (0,0)
    while all(all(abs(xy)<=abs(boundaries)))
        xy = W(xy);
        history = [history;xy];
    end
end
plot(history(:,1),history(:,2),"DisplayName","Brownian Motion")
grid on
hold on
scatter(0,0,'r','filled',SizeData=50)
xlim([-1.2,1.2])
ylim([-1.2,1.2])
xline(boundaries(:,1),'r')
yline(boundaries(:,1),'r')
hold off
%% Task B) Monte Carlo Estimation
dt = 0.0001; % Specified in the assignment.
sigma = 1;
W = @(p) p+sqrt(dt)*sigma*randn(size(p)); % Function to update position.
boundaries = [-1,-1;1,1];
N=10000;
counter = zeros(N,1);
for i=1:N
    xy = zeros(1,2); % Each simulation starts at (0,0)
    while all(all(abs(xy)<=abs(boundaries)))
        xy = W(xy);
        counter(i)=counter(i)+1;
    end
end

%%
mu = mean(counter * dt); % Expected amount of time:
se = std(counter * dt)/sqrt(N); % Standard error: .
% Central Limit Theorem:
% steps~N(mu,se^2) =>
% steps95 = mu+-1.96*se
steps95 = [mu-1.96*se,mu+1.96*se]
% [597.1 - 650.5] dt = 10^-3
% [57310 - 62878] dt = 10^-4
% [57719 - 62918] dt = 10^-5

%% Task C) Monte Carlo Estimation
N = 10000;
dt = 10^-4; % Time Steps
sigs = 1:10;
counter = zeros(N,length(sigs));
k = 0;
for sigma = sigs
    k = k+1;
    W = @(p) p+sqrt(dt)*sigma*randn(size(p));
    for i=1:N
        xy = zeros(1,2); % Each simulation starts at (0,0)
        while all(all(abs(xy)<=abs(boundaries))) % Check if outside any.
            xy = W(xy);
            counter(i,k)=counter(i,k)+1;
        end
    end
end
%%
mu_sigs = mean(counter*dt); % Expected amount of time.
se_sigs = std(counter*dt)/sqrt(N);
% Interesting seems to be mu_sigs/se_sigs = Const.
se_sigs./mu_sigs

plot(sigs, mu_sigs)
xlabel('\sigma')
ylabel('mean')

figure
plot(log10(sigs), log10(mu_sigs))
xlabel('log_{10}(\sigma)')
ylabel('log_{10}(mean)')

