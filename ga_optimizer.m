% Keep codespace clean to avoid bugs and improve readability
clear; clc;
% System
[Phi_system, Phi_q, nu, gamma] = lib.setup_constraints();
nb = 6;

% Options
options = struct();
options.normCorrection = 1.0;
options.tolerance = 1e-8;
options.maxIter = 1000;

% Initial Guess
[Q0, L0] = lib.get_inital_pos(Phi_system, Phi_q, nu, gamma, nb, options);
Lw = 0.85*L0;

% Omega grid
grid_size = 9;
omegas = lib.get_omega_pairs(grid_size);

% Time
time = struct();
time.steps = 100;
time.t0 = 0;

% Setup Genetic algorithm
tic % Time start

% Upper and lower search bounds
L_lb = 0.8*L0;
L_ub = 1.2*L0;
nVars = numel(L0);

% Define function for ga
objfun = @(Lw) objective(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, options);

% ga options
opts = optimoptions('ga', ...
    'PopulationSize', 50, ...
    'MaxGenerations', 1000, ...
    'UseParallel', true, ...
    'Display', 'iter');

% Optimize using generic algorithm
[L_opt, fval] = ga(objfun, nVars, [], [], [], [], L_lb, L_ub, [], opts);
toc % End time

% Display results
disp('Optimal lengths:');
disp(L_opt);
disp('Objective function value:');
disp(fval);

%objfun = objective(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, options)
