s0 = 100;
sigma = 0.3;
r = 0.01;
v = 100;
t = [1,2,3,10];%,20];
mu = (r-sigma.^2/2).*t+log(s0);
variance = sigma.^2.*t;
val = log(v);%/s0);
%val = ones(size(mu))*val;
val = (val-mu)./sqrt(variance);
p = normcdf(val,0,1)
p_exc = 1-p

expected = s0.*exp(r.*t)
tot_var = s0.^2.*exp(2*r.*t).*(exp(sigma.^2.*t)-1)