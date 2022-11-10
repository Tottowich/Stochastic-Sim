close all
clear all
mu = 800;
N = 1000;
p = 0.05;
limit_check = 50000;
lamb= 1/mu;
g = @(x)lamb*exp(-lamb*x);
F_inv = @(u)-(1/lamb)*log(1-u);
iter = 1000;
X = F_inv(rand(iter,N)); % Sample the value of the "gift"
u = rand(N,1);
n = binoinv(u,N,p); % Sample the value of
Tot = n'.*mean(X);
exceeding = zeros(size(Tot));
exceeding(Tot>limit_check) = 1; 
% Logical array will be binomaly distributed
p_exc = mean(exceeding);

% The Interval of a binomally distributed variable can be calculated as:
% p_bar+-z*sqrt(p_bar(1-p_bar)/N).
alpha = 0.05;
z = 1.96; % Retrieved from table from normal distribution.
std_exc = sqrt(p_exc*(1-p_exc)/N);
L = p_exc-z*std_exc;
U = p_exc+z*std_exc;
% Normal distribution with mu=p_exc and std=std_exc
x = [p_exc-5*std_exc:0.0001:p_exc+5*std_exc];
y = normpdf(x,p_exc,std_exc);
figure(2)
plot(x,y,"DisplayName","Normal Approximation")
hold on
bar(x(x>=L&x<=U),y(x>=L&x<=U),"FaceColor","b","EdgeColor","b","DisplayName","95%")
grid on
legend()
hold off





