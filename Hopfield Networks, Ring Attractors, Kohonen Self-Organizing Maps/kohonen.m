%% HW10 - 3. Computer experiment with a Hopfield network

%% 1 - Self-organization
% When the number of trials is set very high, the neurons in the network
% converge to an equilibrium state of a preferred direction to fire in
% given an input. In self-organizing maps, a neuron is made a "winner," in
% which that neuron's weight vector matches most closely to the input
% vector, and the neighborhood of the winner neuron updates their weight
% vectors to the same direction as the winner. With 1000 random stimuli,
% neurons that are in close proximity under a Gaussian Distribution are
% tuned to the winner. Thus, the behavior of each neuron will not change
% as much of the neurons in the neighorhood of the winner will already be
% pointing in the same direction.

%%
SOM(1000, 1.5, 2)
SOM(1000, 1.5, 2)

%% 2 - Narrow neighborhood
% With a narrower Gaussian neighborhood, the winner neuron from each random
% stimuli does not influence its neighbors very much, so as a result, the
% network does not converge to a uniform field that prefers a particular
% input. This makes sense, because less neurons will be in the vicinity of 
% the winner neuron, and thus less neurons are tuned/adjusted to have a 
% specific weight vector at very trial.

%%
SOM(1000, 0.1, 2)
SOM(1000, 0.1, 2)

%% 3 - Wide neighborhood
% With a wider Gaussian neighborhood, the winner neuron in every trial
% influences almost all of the neurons in the network to have the same
% weight vector, which I observed as different inputs were given to the
% network.

%%
SOM(1000, 1000, 2)
SOM(1000, 1000, 2)

%% 4 - Restricted inputs
% By restricting inputs to be only in the top horizontal half of the grid,
% neurons are tuned to have weight vectors that point also to the top 
% horizontal half of the grid. This is because of the neurons that are
% chosen as winners are generally pointing in upwards, thus tuning other
% neurons in its neighborhood to point upwards as well. After a high number
% of inputs in only the upper half grid, neurons pointing in the lower half
% grid will lose their selectivity of inputs pointing downwards.

%%
SOM(1000, 1.5, 1)
SOM(1000, 1.5, 1)
