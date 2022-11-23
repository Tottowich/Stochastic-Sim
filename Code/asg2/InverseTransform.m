function x=InverseTransform(p)
    % Given at Lecture "Markov Chains"
    F=cumsum(p);
    u=rand;
    x=find(u<F,1,'first');
end