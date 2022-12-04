%%
p = 0.5;
steps = 10;
sims = 2000000;
pos = -5;
disp("Simulating Random Walk")
rand_walk = RandomWalk(steps,p,sims);
present=any(rand_walk==pos,1);
mean(present)