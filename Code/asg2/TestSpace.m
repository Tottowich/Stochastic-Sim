%% Testing various functions.
intensity = @(t) 30.*t.*(1-t);
lambda_max = 7.5;
T = 1;
N = 1000000;
sims = zeros(N,1);
for i = 1:N
    sims(i) = length(nonhomogenous(intensity,lambda_max,T));
end
mean(sims) % Checks Out!