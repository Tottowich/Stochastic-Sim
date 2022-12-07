%% Task1: Rat Task (Markov Chain)
P = [0 1/2 0 0 0 1/2 0 0 0;
        1/3 0 1/3 0 1/3 0 0 0 0;
        0 1/2 0 1/2 0 0 0 0 0;
        0 0 1/3 0 1/3 0 0 0 1/3;
        0 1/4 0 1/4 0 1/4 0 1/4 0
        1/3 0 0 0 1/3 0 1/3 0 0;
        0 0 0 0 0 1/2 0 1/2 0;
        0 0 0 0 1/3 0 1/3 0 1/3;
        0 0 0 1/2 0 0 0 1/2 0;
    ]; % The given transition matrix
%% Rat Task A)
counter = zeros(1,9);
N = 100000;
%for k = 1:9
    p = 1; %k;
    for i = 1:N
        pos_vec = zeros(1,9);
        pos_vec(p) = 1;
        counter(p) = counter(p) + 1;
        pos_vec = pos_vec*P;
        p = InverseTransform(pos_vec);
    end
%end
dist = counter/sum(counter)

pos_vec = zeros(1,9);
pos_vec(1) = 1;
dist_1 = pos_vec*P^100000 % Test if starting in odd state

pos_vec = zeros(1,9);
pos_vec(2) = 1;
dist_2 = pos_vec*P^100000 % Test if starting in odd state
%% Rat Task B)

% Calculate Steady state:
% The matrix is irridusabel and aperiodic
% NB - not aperiodic

% You can get from any cell to any
% True

% In k and k+1 steps can you go from
% one cell to itself
% THIS IS FALSE
% Only when k is even can you have gotten back (equally many steps left as
% right, and up as down). When odd you can never have gotten to the same
% square!

% Solving this theoretically:

% pi*P = pi => Asymtotic solution to problem if sum(pi)=1.

% pi*P = pi <=> pi*P - pi=0 <=> pi*(P-I)=zeros(1,N)

% Compute P-I but substitute one column with zeros and one entry in
% zeros(1,N). This is to induce that sum(pi)=1.
% A = P-I:
n = size(P,1);
A = P-eye(n);
i = 1;
A(:,i) = 1;  % <== Substituting column
O_h = zeros(n,1);
O_h(i) = 1; % Create a one hot encoded right hand side.
% Then solve: pi*A= O_h
pi_vec = A'\O_h % Theoretical Distribution


%% Rat task C)
N = 1000;
counter = zeros(N,1);
pos_vec = zeros(1,9);%randi(9);
pos_vec(5) = 1; % Start at pos 5.
for i = 1:N
    p = 5;
    while p~=5 || counter(i)==0
        pos_vec = zeros(1,9);
        pos_vec(p) = 1;
        counter(i) = counter(i) + 1;
        pos_vec = pos_vec*P;
        p = InverseTransform(pos_vec);
    end
end
% According to the central limit theorem it can be normally approximated since N>30 
% with:
%   smu = mean(counter) and sigma = std(counter)/sqrt(N)
mu = mean(counter)
sigma = std(counter)/sqrt(N) % Standard error.


