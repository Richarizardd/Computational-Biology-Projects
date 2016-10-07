%% HW10 - 2. Computer experiment on a ring network

%% 1 - Stationary activity bump
% As seen below, the system does now always reach the same equilibrium
% state. This is because random initial states were chosen, and as a
% result, the time evolution of these states converges to different equilibrium states
% , which represents different sets of neurons being fired.

%% 
% B = -1.5, S = 0.0

%%
ringnet(-1.5, 0, 0)
ringnet(-1.5, 0, 0)

%% 2 - Traveling activity bump
% With S = 0, our ring network has fairly symmetric lateral connections, 
% and thus do not change states. However, by having S be non-zero, we see a 
% "rotation" of neurons firing at different time points in a loop, which is
% due to asymmetric lateral connections. By y increasing the magnitude of 
% the weight shift parameter, our network change states more rapidly, and 
% thus, "oscillate" at a greater frequency, with the sign of S indicating 
% the direction of asymmetry/states are traverse.

%%
% B = -1.5, S = -0.5, No Noise

%%
ringnet(-1.5, -0.5, 0)

%%
% B = -1.5, S = -0.2, No noise

%%
ringnet(-1.5, -0.2, 0)

%%
% B = -1.5, S = 0.0, No noise

%%
ringnet(-1.5, 0.0, 0)

%%
% B = -1.5, S = 0.2, No noise

%%
ringnet(-1.5, 0.2, 0)

%%
% B = -1.5, S = 0.5, No noise

%%
ringnet(-1.5, 0.5, 0)

%% 3 - Unfiform state vs symmetry breaking
% By reducing the inhibitory level of B, the system saturates to a maximum
% neuronal activity value of approximately 10 across all time points. 
% Essentially, our neurons do not have a preference for any direction, as
% all neurons become excited.

%%
% B = -1.0, S = -0.5, No noise

%%
ringnet(-1.0, -0.5, 0)

%%
% B = -1.0, S = -0.2, No noise

%%
ringnet(-1.0, -0.2, 0)

%%
% B = -1.5, S = 0.0, No noise

%%
ringnet(-1.0, 0.0, 0)

%%
% B = -1.0, S = 0.2, No noise

%%
ringnet(-1.0, 0.2, 0)

%%
% B = -1.0, S = 0.5, No noise

%%
ringnet(-1.0, 0.5, 0)

%% 4 - Noisy weights
% In adding noisy weights to the weight matrix, the neurons in our network 
% do not exhibit the same pattern of lateral connections. Thus, there is 
% not a clear pattern of r ecurrent weight distribution. In the case of S =
% 0, we see that our network still converges to an equilbrium state, 
% however, the level of excitation varies drastically. For S not equal to 0
% , roughly, we still see rotation-invarance reflected in the activity, as 
% the peak intensities still shift from row to row, but the rotation is not
% very smooth.

%%
% B = -1.5, S = -0.5, With noise

%%
ringnet(-1.5, -0.5, 1)

%%
% B = -1.5, S = -0.2, With noise

%%
ringnet(-1.5, -0.2, 1)

%%
% B = -1.5, S = 0.0, With noise

%%
ringnet(-1.5, 0.0, 1)

%%
% B = -1.5, S = 0.2, With noise

%%
ringnet(-1.5, 0.2, 1)

%%
% B = -1.5, S = 0.5, With noise

%%
ringnet(-1.5, 0.5, 1)