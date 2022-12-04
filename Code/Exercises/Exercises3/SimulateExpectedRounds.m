
% Program so simulate the expected number of rounds from of a "Random Game"

% a - Starting amount of the person with probability p to win one round.
B = 100;    % Total Balance of the game.
sims = 1000;
steps = 6000;
p = 0.8;
a_s = 1:B-1;
simulations = zeros(size(a_s,2),2);
for a=a_s
    % Total Balance of the game.
    b = B-a;
    % Check how long it takes for a to exceed B or fall below 1.
    rand_walk = RandomWalk(steps,p,sims,a);
    endings = rand_walk==0 | rand_walk>=B;
    [redun,first_zero]=max(endings);
    fprintf("Simulated number of steps for a=%d.\n",a)
    simulated = mean(first_zero);
    [P,expected] = RandomGame(p,a,b,1);
    simulations(a,1)=expected;
    simulations(a,2)=simulated;
end


%%
% Display the maximum expected number of rounds.
[max_expected, max_a] = max(simulations(:,1));
[min_expected, min_a] = min(simulations(:,1));
[min_sim, min_sim_a] = min(simulations(:,2));
[max_sim, max_sim_a] = max(simulations(:,2));

fprintf("Results for %d simulations of p=%.3f and B=%d:\n",sims,p,B)
fprintf("\tMinimum expected number of rounds: %4.2f, for a=%d.\n",min_expected,min_a);
fprintf("\tMaximum expected number of rounds: %4.2f, for a=%d.\n",max_expected,max_a);
fprintf("\n\tMinimum simulated number of rounds: %4.2f, for a=%d.\n",min_sim,min_sim_a);
fprintf("\tMaximum simulated number of rounds: %4.2f, for a=%d.\n",max_sim,max_sim_a);

% Plot the results of the simulation vs the expected number of rounds.

figure(1)
plot(simulations(:,1),'r')
hold on
plot(simulations(:,2),'b')r
hold off
grid on
xlabel('a')
ylabel('Number of rounds')
legend('Expected','Simulated')
title(sprintf('Expected vs Simulated number of rounds for p=%.3f and B=%d',p,B))







