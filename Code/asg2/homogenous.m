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