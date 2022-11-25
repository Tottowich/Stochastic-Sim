function S = nonhomogenous_decomposed(intensity, lambda_max, T)
% NONHOMOGENOUS_DECOMPOSED  Simulate non-homogenous Poisson process events.
%
% CALL SEQUENCE:
%
%   [ S ] = nonhomogenous_decomposed(intensity, lambda_max, T)
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
%   This is the same as the thinning algorithm, but maybe a bit clearer
%   what is being done.
%
% MINIMAL WORKING EXAMPLE:
%
%   T = 10;
%   intensity = @(t) 2 + (3./sqrt(t + 1));
%   lambda_max = 5;
%   S = nonhomogenous_decomposed(intensity, lambda_max, T);
%   N = length(S);
%

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