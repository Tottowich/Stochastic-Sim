%% Task3: Decomposition (Poisson process)

% Consider a non-homogenous Poisson process {N(t), t \gte 0} 
% with intensity function \lambda(t) = 2 + 3/sqrt(t+1), t>0

% We want to simulate this process on the time interval [0,10]

% We have a homogenous Poisson process of events  occuring at
% random time points, with N(t) denoting denoting the number of events
% in interval [0,t), if:
%   - N(0) = 0
%   - The increments are independent, meaning that the number of events 
%       occuring in disjoint intervals are independent.
%   - The increments are stationary, meaning that the distribution of N
%       depends only the length of the interval, and not its location.
% 
% The number of occurances of events of a Poisson process with intensity
% \lambda in a time interval of t > 0 is Poisson distributed with parameter
% t\lambda. That is N(s+t) - N(s) \sim Po(t\lambda), \forall s>0.
% (We have the special case when s = 0 that N(0) = 0 so 
% N(t) \sim Po(t\lambda).)
%
% A Poisson distribution with intensity parameter \lambda has 
% P(X = k) = ((\lambda)^k / k!) e^{-\lambda}
% and has expected value and variance \lambda.
%
% A nonhomogenous Poisson process has intensity varying with time,
% \lambda(t) \gte 0. It can better model events where the dynamics of arrival
% depents on time. (It will have sligthly different properties because of this.)
% 
% To find the distribution of the number of events of a non-homogenous
% Poisson process we begin by defining the mean-value function \Lambda(t):
% \Lambda(t) = \integral_0^t \lambda(\tau)d\tau
% Then the number of occurances in interval (s,s+t) is Poisson distributed
% with parameter \Lambda(s+t) - \Lambda(s):
% N(s+t) - N(s) \sim Po(\Lambda(s+t) - \Lambda(s))
% ( Or in our case N(t) \sim Po(\Lambda(t)-\Lambda(0) ).
%
% One way to create a non-homogenous process with intensity \lambda(t)
% from a homogenous one with intensity \lambda, is to only let through
% an event at time t with probability p(t) = \lambda(t)/\lambda. (Given the
% original process \lambda \gte max_{t \in [0,T]} \lambda(t) .)
%


%% a Decompose and compute homogenous and nonhomogenous components
help("homogenous")
help("nonhomogenous")

%% Homogenous process algorithm sanity check
lmb = 0.42;
T = 100;

N = 10000;

res_ours = zeros(N,1);
res_theirs = zeros(N,1);
for i = 1:N 
    S = homogenous(lmb, T);
    res_ours(i) = length(S);
    
    res_theirs(i) = poissrnd(lmb*T);
end 

fprintf('Sanity checking using lambda=%.2f, T=%d\n', lmb, T)
fprintf('Built-in poissrnd result: %f\n', mean(res_theirs))
fprintf('Simulated homogenous process: %f\n', mean(res_ours))
fprintf('Expected intensity: %.2f\n', T*lmb)

%% Non-Homogenous process thinning algorithm sanity check

T = 10;
intensity = @(t) 2 + (3./sqrt(t + 1));

meanvaluefunc = @(t) 2*t + 6 * sqrt(t + 1);
lambda_tot = meanvaluefunc(T) - meanvaluefunc(0);

% lambda-max for t \in [0,10] for \lambda(t) = 2 + 3/sqrt(1+t)
% is given at t = 0, where \lambda(t) = 5.
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

fprintf('Sanity checking using given intensity function for T=%d\n', T)
fprintf('Built-in poissrnd result: %f\n', mean(res_theirs))
fprintf('Simulated homogenous process: %f\n', mean(res_ours))

%% Decompose Non-Homogenous process

% The intensity of the combination of the events from two Poisson processes
% is is the sum of the two separate process intensities.
% In our case we can treat 2 + 3/sqrt(1+t) as the sum of one homogenous
% process with intensity \lambda_h = 2; and one non-homogenous process with
% intensity \lambda_n(t) = 3/sqrt(1+t).

intensity_h = 2;
intensity_n = @(t) (3./sqrt(t + 1));

lambda_max_n = 3;
meanvaluefunc_n = @(t) 6 * sqrt(t + 1);

simulate = @(t) sort( [ ...
    homogenous(intensity_h, t), ...
    nonhomogenous(intensity_n, lambda_max_n, t) ...
    ] );


%% b Simulate trajectories

max_T = 10;
ts = 0:0.1:max_T;

N = 50;
data_N = zeros(N, length(ts));
data_S = cell(N, length(ts));

means_theory = zeros(length(ts), 1);
for i = 1:N
    K = 0;
    for T = ts 
        K = K + 1;

        S = simulate(T);
        data_S{i, K} = S;

        Nt = length(S);
        data_N(i, K) = Nt;

        if i == 1
            means_theory(K) = meanvaluefunc(T) - meanvaluefunc(0);
        end
    end
end

% Plot trajectories
max_T_ind = length(ts);

figure
hold on
for i = 1:N
    Nt = data_N(i,max_T_ind);
    S = data_S{i,max_T_ind};
    
    stairs([0, S, max_T], [0:Nt, Nt], ...
        LineWidth=0.1, HandleVisibility="off")
end

plot(ts, meanvaluefunc(ts)-meanvaluefunc(0), ...
    'k', LineWidth=3, DisplayName='\Lambda(t)-\Lambda(0)')
legend(Location="best")
xlabel('T')
ylabel('N(T)')
title('Stair plot for nonhomogenous poisson process')
hold off

% Plot number of events
means_approx = mean(data_N, 1);
figure
hold on
grid on

plot(ts, data_N, ...
    HandleVisibility="off", LineWidth=0.1, ...
    ColorMode="auto", Color=[rand(3,1); 0.25])

plot(ts, means_theory, ...
    LineWidth=3, DisplayName='means')
plot(ts, means_approx, ...
    LineWidth=3, DisplayName='approx')

legend(Location="best")
xlabel('T')
ylabel('N(T)')
title('Number of events of nonhomogenous poisson process')
hold off

%% c Estimate expected value and variance of N(10)

N = 100000;
data_N = zeros(N, 1);
T = 10;
for i = 1:N
    Nt = length(simulate(T));
    data_N(i) = Nt;
end

fprintf('Approximative mean: %f\n', mean(data_N))
fprintf('Approximate variance: %f\n', var(data_N))

fprintf('Theoretical variance and mean: %f\n', meanvaluefunc(T) - meanvaluefunc(0))


%% Homogenous Poisson Process Simulation
% Following algorithm in  Ross, S M., Simulation and course slides.

function S = homogenous(lambda, T)
% HOMOGENOUS  Simulate homogenous Poisson process events.
%
% CALL SEQUENCE:
%
%   [ S ] = homogenous(lambda, T)
%
% INPUT:
%
%   lambda  The intensity of the process. A number > 0.
%   T       The length of time to simulate the process for. A number >= 0.
%
% OUTPUT:
%
%   S       A series of event times, in order.
%
% DESCRIPTION:
%
%   Simulate a homogenous Poisson process with intensity lambda running for
%   T length of time. Each event time is then added to the vector S.
%   To get the number of events N(T) = length(S) can be used.
%
% REMARKS:
%
%   Following algorithm in  Ross, S M., Simulation and course slides.
%
% MINIMAL WORKING EXAMPLE:
%
%   lmb = 0.42;
%   T = 100;
%   S = homogenous(lmb, T);
%   N = length(S); % Should be approximatively T * lmb.
%

    assert(lambda > 0, "Expect intensity lambda positive")
    assert(T >= 0, "Expect time length positive")

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

%% Nonhomogenous Poisson Process Simulation - Thinning algorithm
% Following algorithm in  Ross, S M., Simulation and course slides.

function S = nonhomogenous(intensity, lambda_max, T)
% NONHOMOGENOUS  Simulate non-homogenous Poisson process events.
%
% CALL SEQUENCE:
%
%   [ S ] = nonhomogenous(intensity, lambda_max, T)
%
% INPUT:
%
%   intensity   The time-varying intensity function \lambda(t). Must be
%               taking values in the range [0, lambda_max] over the domain [0,T].
%   lambda_max  The intensity of the process. A number > 0, such that the
%               intensity function is bounded by it on the range [0,T].
%   T           The length of time to simulate the process for. A number >= 0.
%
% OUTPUT:
%
%   S           A series of event times, in order.
%
% DESCRIPTION:
%
%   Simulate a non-homogenous Poisson process with intensity-function as
%   given, running for T length of time. Each event time is then added to 
%   the vector S.
%   To get the number of events N(T) = length(S) can be used.
%
% REMARKS:
%
%   Following algorithm in  Ross, S M., Simulation and course slides.
%
% MINIMAL WORKING EXAMPLE:
%
%   T = 10;
%   intensity = @(t) 2 + (3./sqrt(t + 1));
%   lambda_max = 5;
%   S = nonhomogenous(intensity, lambda_max, T);
%   N = length(S);
%

    assert(lambda_max > 0, "Expect lambda_max positive")
    assert(T >= 0, "Expect time length positive")

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

%% Nonhomogenous Poisson Process Simulation - Decomposition algorithm
% This is the same as the thinning algorithm, but maybe a bit clearer
% what is being done.

function S = nonhomogenous_decomposed(intensity, lambda_max, T)

    % given a homogenous process with parameter lambda_max
    S_homogenous = homogenous(lambda_max, T);

    S = zeros(length(S_homogenous), 1);
    I = 0;

    for i = 1:length(S_homogenous)
        t = S_homogenous(i);

        % include event with correct proportion lambda(t)/lambda_max
        if rand <= intensity(t) / lambda_max
            I = I + 1;
            S(I) = t;
        end
    end

    S = S(1:I);
end