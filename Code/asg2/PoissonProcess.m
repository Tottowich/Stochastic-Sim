%% Task3: Decomposition (Poisson process)

% poisson distribution  P(X = k) = (\lambda^k / k!) e^{-\lambda}


%% a Decompose and compute homogenous and nonhomogenous components


% Homogenous

lmb = 5
T = 100

N = 10000

res_ours = zeros(N,1);
res_theirs = zeros(N,1);
for i = 1:N 
    S = homogenous(lmb, T);
    res_ours(i) = length(S);
    
    res_theirs(i) = poissrnd(lmb*T);

end 

mean(res_theirs)
mean(res_ours)

% Non-Homogenous

T = 10;

intensity = @(t) 2 + (3./sqrt(t + 1));

integral = @(t) 2*t + 6 * sqrt(t + 1);
lambda_tot = integral(T) - integral(0);

% lambda-max for t \in [0,10]
% ...
lambda_max = 5;

S = nonhomogenous(intensity, lambda_max, T);
length(S);
poissrnd(lambda_tot);

res_ours = zeros(N,1);
res_theirs = zeros(N,1);
for i = 1:N 
    S = nonhomogenous(intensity, lambda_max, T);
    res_ours(i) = length(S);
    
    res_theirs(i) = poissrnd(lambda_tot);

end 

mean(res_theirs)
mean(res_ours)

%% b Simulate trajectories

ts = [0:0.1:10];

N = 50;
data = zeros(N, length(ts));
means_theory = zeros(length(ts), 1);
for i = 1:N
    K = 0;
    for T = ts 
        K = K + 1;
        Nt = length(nonhomogenous(intensity, lambda_max, T));
        data(i, K) = Nt;

        if i == 1
            means_theory(K) = integral(T) - integral(0);
        end
    end
end

means_approx = mean(data, 1);

plot(ts, data, HandleVisibility="off", LineWidth=0.1, ColorMode="auto", Color=[rand(3,1); 0.25])

hold on
grid on

plot(ts, means_theory, LineWidth=3, DisplayName='means')
plot(ts, means_approx, LineWidth=3, DisplayName='approx')
legend(Location="best")
xlabel('T')
ylabel('N(T)')
title('Number of events of nonhomogenous poisson process')
hold off

%% c Estimate expected value and variance of N(10)

N = 100000;
data = zeros(N, 1);
T = 10;
for i = 1:N
    Nt = length(nonhomogenous(intensity, lambda_max, T));
    data(i) = Nt;
end


plot(1:N, data, HandleVisibility="off")

grid on

hold off

fprintf('Approximative mean: %f\n', mean(data))
fprintf('Approximate variance: %f\n', var(data))

fprintf('Theoretical variance and mean: %f\n', integral(T) - integral(0))


%% Homogenous

function S = homogenous(lambda, T)

    t = 0;
    I = 0;
    S = [];

    while true
        I = I + 1;

        U = rand;
        X = -(1/lambda) * log(U);

        t = t+X;

        if t > T
            break
        end

        S(I) = t;
    end
end

function S = nonhomogenous(intensity, lambda_max, T)

    t = 0;
    I = 0;
    S = [];

    while true

        U = rand;

        t = t - ((1/lambda_max) * log(U));

        if t > T
            break
        end

        U = rand;
        if U <= intensity(t) / lambda_max
            I = I + 1;
            S(I) = t;
        end
    end

end