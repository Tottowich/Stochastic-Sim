% Task 4.
e = exp(1);
f = @(x) exp(x)/(e-1);
%u = rand;
F_inv = @(u) log((e-1)*u+1);
u_s = [0.76,0.74,0.39,0.66,0.17];
x_s = F_inv(u_s);
%%
% Task 4 + acceptance rejection.
g = @(x) 1;
G_inv = @(u) u;
% Determine c as the maximum of max(f/g)=max(f)=e/(e-1)
c = e/(e-1);
N = 10000;
X = zeros(1,N);
for i = 1:N
    while true
        u = rand;
        y = G_inv(u);
        u = rand;
        if c*g(y)*u<f(y)
            x=y;
            break
        end
    end
    X(i)=x;

end
histogram(X,"Normalization","pdf","BinWidth",0.05)
hold on;
t = [0:0.01:1];
plot(t,f(t),'LineWidth',3)
hold off
%%
% Acceptance rejection is the most simple since it has 
f = @(x) 30*(x.^2-2*x.^3+x.^4);
% Derive the function and set the derivative = 0.
% Solve for x and calculate the maximum.
% Max at x = 0.5
% Since only generating from 0->1 => Let:
g = @(x) 1;
% c = max(f/g) => c = f(0.5) = 30*(1/4-1/4+1/16) = 30/16
c = 30/16;
% Inverse cumnatative function =>
G_inv = @(u) u;

% A-R method
N = 100000;
X = zeros(1,N);
for i=1:N
    while true
        u = rand; % Uniform 0->1
        y = G_inv(u);
        u = rand;
        if c*g(y)*u<f(y)
            x=y;
            break
        end
    end
    X(i)=x;
end
t = linspace(0,1,N);
plot(t,f(t))% Probability mass function to emulate.
hold on
grid on
histogram(X,"Normalization","pdf")

% Correct!

%%

% Generate numbers from weibull distribution wiht delta=2 beta=3

% Formula: f(x) = beta/delta*(x/delta)^(k-1)*e^-
delta = 2;
beta = 3;
f = @(x) beta/delta*(x/delta).^(beta-1).*exp(-(x/delta).^beta);
F = @(x) 1-exp(-(x/delta).^beta);
% Solving for x from F we can use the inverse transform method.
F_inv =@(u) delta*(-log(1-u)).^(1/beta);
N = 1000;
X = F_inv(rand(1,N));
t = linspace(1,100,N);
plot(t,f(t))% Probability mass function to emulate.
hold on
grid on
histogram(X,"Normalization","pdf")
hold off;
% Get probabilty that X<2

Y = zeros(size(X));
Y(X<2)=1;
mean(Y)
%%
clear all
close all
possible = 2:12;
reps = 1000;
rolls = zeros(1,reps);
die = @() floor(rand*6+1);
disp("Starting")
founds = zeros(size(possible))
for rep = 1:reps
    found = zeros(size(possible));
    %while prod(ismember(found,possible))
    r = 0;
    while ~prod(found)
        r = r+1;
        d1 = die();
        d2 = die();
        s = d1+d2-1;
        found(s)=found(s)+1;
    end
    rolls(rep) = r;
    founds = founds + found;
end
figure(1)
histogram(rolls,"Normalization","pdf")
mean(rolls)

histogram(founds,"Normalization","pdf")

%%




