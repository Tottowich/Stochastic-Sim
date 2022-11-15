
% Task 1 b.
% Simulate from negative binomial distribution using recurrence.
% We use recurrence since it is r in eq (1.b.1)
N = 1000;
X = zeros(1,N);
R = 10;
P = 0.5
p_next = @(p_j,j) j.*(1-P)./(j+1-R).*p_j;
for i = 1:N
    u = rand;
    j = R;
    p_j = P^R;
    F = p_j;
    while F<u
        p_j = p_next(p_j,j);
        F = F+p_j;
        j = j+1;
    end
    X(i) = j; %
end

