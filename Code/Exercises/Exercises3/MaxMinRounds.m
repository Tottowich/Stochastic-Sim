function a = MaxMinRounds(p,B)
    % MaxMinRounds: Find the maximum number of rounds for a given probability and total balance of a "Random Game"
    % Input:
    %   p: probability of winning for player with balance a
    %   B: total balance of the game
    % Output:
    %   a: Vector containing the a giving the maximum and minimum number of rounds. [max a, min a].
    %ab = a*(B-a);
    %if p == 1/2
    %    a = [ab,ab];
    %else
    Q = (1-p)/p;
    a = log(Q^B-1)-log(abs(log(Q)))-log(B);
end
        
