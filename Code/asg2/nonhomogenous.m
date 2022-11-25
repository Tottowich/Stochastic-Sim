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