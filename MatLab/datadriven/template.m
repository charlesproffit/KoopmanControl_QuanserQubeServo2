%% Template File

% Controller order
order = ...

% Number of measurements
ny = ... 

% Number of actuators
nu = ...

% Initial stabilising controller
K_init = ...

% Frequency Grid
omegas = ...

%% LFR Models
G = ...
G = mktito(G, ny, nu);

%% Problem Formulation

% Controller Object 
rng(22);
K = datadriven.Controller.SS(order, ny, nu, Ts);
K.setinitial(K_init);

% Synthesiser Object
synthesiser = Synthesiser(K);

% Add Hinf objectives
% synthesiser.add_Hinf_objective(G, omegas);

% Add H2 objectives
% synthesiser.add_H2_objective(G, omegas);

% Add Hinf constraints
% synthesiser.add_Hinf_constraint(G, omegas);

% Add H2 objectives
% synthesiser.add_H2_constraint(G, omegas);

% Ensure that controller is stable
synthesiser.ensure_controller_stability(omegas);

%% Optimise
output = synthesiser.synthesise();

K_final = output.Controller;