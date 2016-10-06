% Authors: Scott Albert, Andrew Long
% Learning Theory Final Project 2015
% Summary: This script is the outermost script that you will use to
%          simulate all 5 cases described in the problem statement.

%% Parameter specification

% Time constants
tau1eye = 0.224;    % seconds
tau2eye = 0.013;    % seconds
alpha1eye = 0.004;  % seconds
tau1head = 0.270;   % seconds
tau2head = 0.15;    % seconds
alpha1head = 0.010; % seconds

% Covariances matrices for state and motor updates
Qx = 10^(-5) * eye(7);
Qy = 0.01 * eye(2);

% The goal direction
goal = 40 * pi/ 180;

% The H matrix
H = [-1,0,0,-1,0,0,1; 1,0,0,0,0,0,0];

% Specify the length of each time step
Delta = 0.002;  % seconds

% Specify the total simulation period
p = 0.3; % seconds

% The time at which we want our gaze to are at the target
k1 = 0.1; % seconds

% The number of discrete steps to reach time point k1
k1Steps = k1 / Delta;

% The number of time steps in the simulation
pSteps = p / Delta;

% The initial state vector
x0 = [0;0;0;0;0;0;goal]; 

% Specify the value of the T matrix. Note that T has this value only after
% k1 time has elapsed. It is zero before.
T =  [10^5,0; 0,300];

% Specify the matrix that determines the motor commands cost
L = [10,0;0,10];

% Number of simulations to average
numAverage = 50;

%% Generate A and B matrices
[A, B] = calculateAandB(tau1eye, tau2eye, alpha1eye, tau1head,tau2head, alpha1head , Delta);
disp(A); disp(B);

%% Case 1: No Signal Dependent Noise, No Hold
% Uncertainty Matrices
C1 = [0,0; 0,0];
C2 = [0,0;0,0];

% The length of the holding period
Thold = 0;

% Perform the simulations and plot the results
[K1,G1] = simulateCase(x0, A,B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, goal, 1, 2, ...
    'Case 1: No Signal Dependent Noise, No Hold');

%% Case 2: Small Signal Dependent Noise, No Hold
% Uncertainty Matrices
C1 = [0.01,0; 0,0];
C2 = [0,0;0,0.01];

% The length of the holding period
Thold = 0;

% Perform the simulations and plot the results
[K2,G2] = simulateCase(x0, A,B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, goal, 3, 4, ...
    'Case 2: Small Signal Dependent Noise, No Hold');

%% Case 3: Large Signal Dependent Noise, No Hold
% Uncertainty Matrices
C1 = [2.00,0; 0,0];
C2 = [0,0;0,2.00];

% The length of the holding period
Thold = 0;

% Perform the simulations and plot the results
[K3,G3] = simulateCase(x0, A,B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, goal, 5, 6, ...
    'Case 3: Large Signal Dependent Noise, No Hold');

%% Case 4: Small Signal Dependent Noise, Short Head Hold
% Uncertainty Matrices
C1 = [0.01,0; 0,0];
C2 = [0,0;0,0.01];

% The length of the holding period
Thold = 0.050 / Delta;

% Perform the simulations and plot the results
[K4,G4] = simulateCase(x0, A,B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, goal, 7, 8, ...
    'Case 4: Small Signal Dependent Noise, Short Head Hold');

%% Case 5: Small Signal Dependent Noise, Long Head Hold
% Uncertainty Matrices
C1 = [0.01,0; 0,0];
C2 = [0,0;0,0.01];

% The length of the holding period
Thold = 0.100 / Delta;

% Perform the simulations and plot the results
[K5,G5] = simulateCase(x0, A,B, pSteps, k1Steps, T, Thold, numAverage, ...
    Qx, Qy, H, C1, C2, L, Delta, goal, 9, 10, ...
    'Case 5: Small Signal Dependent Noise, Long Head Hold');

%% Amalgamating kalman gains and feedback gains
kalmanGains = [K1' K2' K3' K4' K5'];
feedbackGains = [G1' G2' G3' G4' G5']; 
save('gains.mat');