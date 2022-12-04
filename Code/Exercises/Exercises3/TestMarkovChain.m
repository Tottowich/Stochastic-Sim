function pivec = TestMarkovChain(P,pivec)
    arguments
        P (:,:) double {mustBeNonnegative, mustBeFinite, mustBeReal}
        pivec (1,:) double {mustBeNonnegative, mustBeFinite, mustBeReal} = 0
    end
    % Input arguments:
    % P: transition matrix, must be provided.
    % pivec: initial distribution, default is uniform.
    % Output arguments:
    % pivec: Asymtotic distribution
    if pivec == 0
        pivec = rand(1,size(P,1));
        pivec = pivec/sum(pivec);
    end
    % Check that P is a stochastic matrix
    if ~all(sum(P,2) == 1)
        error('P is not a stochastic matrix')
    end
    % Check that pivec is a probability vector
    if ~all(sum(pivec,2) == 1)
        error('pivec is not a probability vector')
    end
    % Check that the dimensions of P and pivec are compatible
    if size(P,1) ~= size(pivec,2)
        error('The dimensions of P and pivec are not compatible')
    end
    
    % Compute the asymptotic distribution
    pivec = pivec*P;
    step_size = 100;
    iter = 1;
    while max(abs(pivec - pivec*P)) > 1e-6 && iter < 10000
        pivec = pivec*P^step_size;
        iter = iter + step_size;
    end
    if iter == 10000
        warning('The algorithm did not converge')
    end
end