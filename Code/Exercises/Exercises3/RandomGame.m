function [P,E] = RandomGame(p,a,b,d)
    % RandomGame: Random game with two players
    % Input:
    %   p: Probability of player 1 winning
    %   a: Starting point for player 1
    %   b: Starting point for player 2
    %   d: Value to transfer from player 1 to player 2 at each step
    % Output:
    %   p: Probability of player 1 winning.
    %   E: Expected value of number of games left.
    a = a/d;
    b = b/d;
    if p==1/2
        E = a*b;
        P = a/(a+b);
    else
        q = 1-p;
        quotient = q/p;
        P = ((quotient^a)-1)/(quotient^(a+b)-1);
        E = a/(q-p)-(a+b)/(q-p)*P;
    end
end