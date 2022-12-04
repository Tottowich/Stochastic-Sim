function pos_vec=RandomWalk(steps,rate,sims,start)
    arguments
        steps int32 {mustBeNonnegative, mustBeFinite, mustBeReal}
        rate double {mustBeNonnegative, mustBeFinite, mustBeReal}
        sims double {mustBeNonnegative, mustBeFinite, mustBeReal} = 1;
        start double {mustBeFinite, mustBeReal} = 0;
    end
    step_vec = rand(steps,sims);%,length(start));
    pos_vec = zeros(steps,sims);%,length(start));
    pos_vec(step_vec < rate)= 1;
    pos_vec(step_vec >= rate) = -1;
    pos_vec = cumsum(pos_vec)+start;
end