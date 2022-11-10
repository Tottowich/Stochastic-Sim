
p = 0.1;
r = 40;
tot = 0;
N = 1000;
X = zeros(1,N);
for i = 1:N    
    tot = 0;
    for j = 1:r
        u = rand;
        x = floor(log(u)/log(1-p))+1;
        tot = tot+x;
    end
    X(i) = tot;
end
histogram(X)
mean(X) % This should be r/p
var(X) % n*p*q, q=1-p
