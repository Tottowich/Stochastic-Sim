M = [1/4 1/2 1/4;1/8 1/2 3/8;1/16 1/4 11/16];
P = [0,0,1];
%s = P*(M^1000);

A = M-eye(size(M));
A(:,1)=1;
R = zeros(1,size(M,2));
R(1)=1;
format rat
pivec = R*inv(A)
format short
%%
M=[0 1/2 1/2; 1/2 0 1/2;1/2 1/2 0]
P=rand(1,3);
P=P/sum(P)
%P=[1/3 1/3 1/3]

P*(M^1000)
%%
M=[0 1/2 1/2; 1/2 0 1/2;1/2 1/2 0]
P=rand(1,3);
P=P/sum(P)
%P=[1/3 1/3 1/3]

P*(M^1000)
%%
% Square transition matrix
clear all
P = [0 1/2 0 1/2; 1/2 0 1/2 0; 0 1/2 0 1/2; 1/2 0 1/2 0];
TestMarkovChain(P,[1/4 1/4 1/4 1/4])

%%

N = 1000000000;
u = rand(1,N);
g = @(u) (1./(u.^2)).*exp(-1./(2*u.^2))./sqrt(2*pi);

mu = mean(g(u));
act = normcdf(-1,0,1);

err = abs(mu-act)

rel_err = err/act

%%

P = [0 0.5 0.5;0.5 0 0.5; 0.5 0.25 0.25]
P^100000


%% Random Game Test


a=(0:10)+3;
b=(10:-1:0)+3;
p = 0.55;
d=1;
[P,E]=RandomGame(p,a,b,d)

%%

a0 = 3;
b0 = 3;
T = 10;
B = a0+b0+T;
q=0.45;
p=0.55;
r = q/p;

c = (log(((r^B)-1)/(B*log(r)-1))/log(r))-a0




