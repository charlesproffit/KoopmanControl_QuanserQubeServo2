clc;
clear;
close all;
addpath(genpath('src'));
addpath(genpath('data'));
addpath(genpath('datadriven'));

rng(43);

% States q :
%   q(1) => Angle of rotary arm
%   q(2) => Angle of pendulum arm
%   q(3) => Angular speed of rotary arm q(1)_dot
%   q(4) => Angular speed of pendulum arm q(2)_dot

%% ------------- SYSTEM PARAMETERS ---------------------------- %%
Ts = 0.008;                                     % Time step
nstates = 4;                                    % Number of states
ninputs = 1;                                    % Number of inputs   
nregions = 16;                                  % Number of subdivisions of a full circle for the initial conditions
angle_region = 2*pi/nregions;                   % Angle range per region for initial angle
ncontrolled_states = 2;                         % For partial FL and collocated control
ncross_val_groups = 5;                          % For ridge regression and choice of lambda
numax = 7;                                      % Max possible input
nsteps = 500;                                   % Total number of steps per trajectory (u(k)=0 if k>n.steps_excited)
ntrajs_training = 160;                          % Trajectories for training EDMD models
ntrajs_testing = 40;                            % Trajectories for testing and simulating control

%% -------- DYNAMICS DEFINITION --------------------------- %%

[Mq,H,Phi,B] = QuanserQubeDynamicsCompactForm();
func_discard = @(x,x0) abs(x(1)) >= pi/2 || abs(x(2)-x0(2)) >= 2*pi || abs(x(2)) >= pi;
Model_NL = NLModel(nstates,ninputs,Ts,Mq,H,Phi,B);
K_LQR = Model_NL.designLQR([0;0;0;0], 0, diag([0.41 0 0 0]), 0.016);

func_lifting = @(q)[
    1;
    q;
    Model_NL.f_nr(q) - Model_NL.g_nr(q)*K_LQR*q;
    Model_NL.g_nr(q);
];
func_linear = @(q)(q);
nlifted_states = size(func_lifting(zeros(nstates,1)),1);

%% --------------DATA COLLECTION----------------------- %%
Sim_NL = OLSimulator(Model_NL,func_discard);

func_initialStates = @() [
    0;
    -pi + angle_region * (randi([2,nregions-1]) + rand() - 1);
    5*pi*(2*rand()-1);
    5*pi*(2*rand()-1)
];
func_inputs = @() 2*numax*rand(ninputs,nsteps) - numax;

% dataset = Sim_NL.generateDataset(ntrajs_training,ntrajs_testing,nsteps,func_initialStates,func_inputs,'euler');
% % OR
[dataset,nsteps,ntrajs_training,ntrajs_testing] = structure_data('data\data.tdms', nstates, ninputs, 0.8, 50, true);


%% -------------- EDMD ----------------------- %%

Model_BILIN = EDMDModel("BILINEAR", "LS", func_lifting, nstates, ninputs, nlifted_states, Ts, ncross_val_groups);
Model_BILIN.EDMDIdentification(dataset);
Sim_BILIN = OLSimulator(Model_BILIN);
X_BILIN = Sim_BILIN.simulate_multistep(dataset.testing.X(:,1,:), dataset.testing.U,'euler');
X_BILIN_onestep = Sim_BILIN.simulate_onestep(dataset.testing.X, dataset.testing.U,'euler');

Model_EDMD = EDMDModel("LINEAR", "LS", func_lifting, nstates, ninputs, nlifted_states, Ts, ncross_val_groups);
Model_EDMD.EDMDIdentification(dataset);
Sim_EDMD = OLSimulator(Model_EDMD);
X_EDMD = Sim_EDMD.simulate_multistep(dataset.testing.X(:,1,:), dataset.testing.U,'euler');
X_EDMD_onestep = Sim_BILIN.simulate_onestep(dataset.testing.X, dataset.testing.U,'euler');

Model_LIN = EDMDModel("LINEAR", "LS", func_linear, nstates, ninputs, nstates, Ts, ncross_val_groups);
Model_LIN.EDMDIdentification(dataset);
Sim_LIN = OLSimulator(Model_LIN);
X_LIN = Sim_LIN.simulate_multistep(dataset.testing.X(:,1,:), dataset.testing.U,'euler');
X_LIN_onestep = Sim_LIN.simulate_onestep(dataset.testing.X, dataset.testing.U,'euler');

data = dataset;
comparison_Identification = compare_EDMD(dataset, X_BILIN, X_EDMD, X_LIN, X_BILIN_onestep, X_EDMD_onestep, X_LIN_onestep, Model_BILIN, Model_EDMD, Model_LIN);

% plots_EDMD(Ts, comparison_Identification);

%% ---------DATA-DRIVEN FEEDBACK LINEARIZATION SYSTEM-------------------- %%

% Reference to track 
q_ref = zeros(nstates,nsteps,ntrajs_testing);

% DDFL
Model_DDFL = DDFLModel(func_lifting,nlifted_states,ncontrolled_states,Ts);
Model_DDFL.DDFLIdentification(Model_BILIN);
Controller_DDFL = FLController(Model_DDFL);
K_DDFL = Controller_DDFL.designLQR(diag([10000, 1]), 1);
Sim_DDFL = CLSimulator(Model_NL,Controller_DDFL);
X_DDFL = Sim_DDFL.generateTrajs(ntrajs_testing,nsteps,func_initialStates,q_ref,'euler');

% MBFL
Model_MBFL = MBFLModel(K_LQR,Ts);
Controller_MBFL = FLController(Model_MBFL);
K_DMBFL = Controller_MBFL.designLQR(diag([10000, 1]), 1);
Sim_MBFL = CLSimulator(Model_NL,Controller_MBFL);
X_MBFL = Sim_MBFL.generateTrajs(ntrajs_testing,nsteps,func_initialStates,q_ref,'euler');

% Aggressive LQR
Controller_AggressiveLQR = SFController(Model_NL);
K_AggressiveLQR = Controller_AggressiveLQR.designLQR([0;0;0;0], 0, diag([10000, 0, 1, 0]), 1);
Sim_AggressiveLQR = CLSimulator(Model_NL, Controller_AggressiveLQR);
X_AggressiveLQR = Sim_AggressiveLQR.generateTrajs(ntrajs_testing,nsteps,func_initialStates,q_ref,'euler');

comparison_Control = compare_FL( {X_DDFL, X_MBFL, X_AggressiveLQR}, {'DDFL', 'MBFL', 'Aggressive LQR'}, Ts, 1000);

%% -------------- RESULTS OF REAL SYSTEM ----------------- %%

[data_results_mbfl, ~,~,~] = structure_data("data\6june\data_MBFL.tdms", nstates, ninputs, 0.8, 500, false);
[data_results_ddfl, ~,~,~] = structure_data('data\6june\data_DDFL.tdms', nstates, ninputs, 0.8, 500,false);
[data_results_lqr, ~,~,~] = structure_data("data\6june\data_lqronly.tdms", nstates, ninputs, 0.8, 500,false);
[data_results_mbfl_track, ~,~,~] = structure_data("data\6june\data_mbfltracking.tdms", nstates, ninputs, 0.8, 'all',false);
[data_results_ddfl_track, ~,~,~] = structure_data('data\6june\data_ddfltracking.tdms', nstates, ninputs, 0.8, 'all',false);
[data_results_agressivelqr_track, ~,~,~] = structure_data("data\6june\data_aggressivelqrtracking.tdms", nstates, ninputs, 0.8, 'all',false);
[data_results_moderatelqr_track, ~,~,~] = structure_data("data\6june\data_moderatelqrtracking.tdms", nstates, ninputs, 0.8, 'all',false);

comparison_Results = compare_FL( ...
    {data_results_lqr.testing, data_results_mbfl.testing, data_results_ddfl.testing}, ...
    {'Aggressive LQR', 'MBFL', 'DDFL'}, Ts, 1000);

comparison_Tracking = compare_FL( ...
    {data_results_mbfl_track.testing, data_results_ddfl_track.testing}, ...
    {'MBFL Tracking', 'DDFL Tracking'}, Ts, 1000);


%% -------------------- ROBUST CONTROL DESIGN --------------------------%%

[data_robust, nsteps,ntrajs_training,ntrajs_testing] = structure_data("data\7june\data.tdms", nstates, ninputs, 1, 'all', true);

%%
% Create controllers object
K_ROB_MIXSYN = RobustController(Model_DDFL);
K_ROB_MIXSYN.buildMultimodelSet(data_robust, 10, func_linear, ncontrolled_states, ninputs, Ts, ncross_val_groups);
%%

% Robustness filter W2
Parray = frd(K_ROB_MIXSYN.G_Set,logspace(-1,2.55,60));
[~,Info] = ucover(Parray,K_ROB_MIXSYN.G_nominal,4);
W2 = Info.W1;

% Performance filter W1
s = tf('s');
W1 = c2d((s+0.0001)/1/(s+0.3),Ts);

% Input filter W3
W3 = c2d(tf(1/18),Ts);
% W3 = [];

figure;
bodemag(K_ROB_MIXSYN.relerr,'b--',W1, 'r', W2, 'g', {0.1, pi/Ts});
legend('Relerr','$W1(s)$','$W2(s)$','Interpreter','latex', 'Location', 'southeast');
title('TFs');
grid on;

%%

% Design mixsyn controller
K_ROB_MIXSYN.designMixsyn(W1, W2, W3);
K_ROB_MIXSYN.validate();
K_ROB_MIXSYN.plotTFs();

%%

K_ROB_DD = RobustController(Model_DDFL);
K_ROB_DD.buildMultimodelSet(data_robust, 10, func_linear, ncontrolled_states, ninputs, Ts, ncross_val_groups);

% Design data-driven robust contoller
K_ROB_DD.designDataDriven(W1, W2, W3, Ts, ncontrolled_states, ninputs);
K_ROB_DD.validate();
K_ROB_DD.plotTFs();

%% 
G_Set = [];
A_sum = zeros(size(Model_DDFL.A_c));
B_sum = zeros(size(Model_DDFL.B_c));

nfolds = 5;
traj_per_fold = floor(ntrajs_training / nfolds);

% Multimodel set
for i = 1:nfolds
    Model_ROBUST = EDMDModel("LINEAR", "LS", func_linear, ncontrolled_states, ninputs, ncontrolled_states, Ts, ncross_val_groups);
    dataset_fold.training.X = data_robust.training.X([1,3], :, traj_per_fold*(i-1)+1 : traj_per_fold*i);
    dataset_fold.training.U = data_robust.training.U(:, :, traj_per_fold*(i-1)+1 : traj_per_fold*i);
    Model_ROBUST.EDMDIdentification(dataset_fold);

    A_sum = A_sum + Model_ROBUST.A_CT;
    B_sum = B_sum + Model_ROBUST.B_CT;

    if isempty(G_Set)
        G_Set = c2d(ss(Model_ROBUST.A_CT, Model_ROBUST.B_CT, [1, 0], 0),Ts);
    else
        G_Set = stack(1, G_Set, c2d(ss(Model_ROBUST.A_CT, Model_ROBUST.B_CT, [1, 0], 0),Ts));
    end
end

% Nominal model (mean of models)
G_nominal = c2d(ss(A_sum / nfolds, B_sum / nfolds, [1, 0], 0),Ts);

% Robustness filter W2
Parray = frd(G_Set,logspace(-1,2.55,60));
relerr = G_nominal\(G_nominal-G_Set);
[P,Info] = ucover(Parray,G_nominal,2);
W2 = Info.W1;

% Performance filter W1
% W1 = makeweight(2,  6, 0.01, Ts);
s = tf('s');
W1 = c2d((s+10)/2/(s+0.000001),Ts);

% Input filter W3
W3 = c2d(tf(1/18),Ts);

figure;
bodemag(relerr,'b--',W1, 'r', W2, 'g', {0.1, pi/Ts});
legend('Relerr','$W1(s)$','$W2(s)$','Interpreter','latex', 'Location', 'southeast');
title('TFs');
grid on;

%% ----------------------Mixsyn------------------------------------------%%
[K_mixsyn,~,gamma] = mixsyn(G_nominal,W1,W3,W2);

% Check if infinity norm if less than 6 dB
infnorm = norm(1/W1, 'inf')
infnorm_dB = 20 * log(infnorm)/log(10)

S = feedback(1, G_nominal * K_mixsyn);
T = feedback(G_nominal * K_mixsyn, 1);
U = feedback(K_mixsyn, G_nominal);

figure;
step(T);
title('Step Response - Output (T) with W3');
grid on;

figure;
step(U);
title('Step Response - Control Signal (U) with W3');
grid on;

figure;
bodemag(S,1/W1);
legend('Sensitivity S(z)','$W1(z)^{-1}$','Interpreter','latex', 'Location', 'southeast');
title('Bode magnitude of sensitivity TF with W3');
grid on;

figure;
bodemag(T,1/W2);
legend('T(z)','$W2(z)^{-1}$','Interpreter','latex', 'Location', 'southeast');
title('Bode magnitude of T with W3');
grid on;

figure;
bodemag(U,1/W3);
legend('Input sensitivity U(z)','$W3(z)^{-1}$','Interpreter','latex', 'Location', 'southeast');
title('Bode magnitude of input sensitivity with W3');
grid on;

[K_mixsyn_TF_b,K_mixsyn_TF_a] = ss2tf(K_mixsyn.A, K_mixsyn.B, K_mixsyn.C, K_mixsyn.D);
K_mixsyn_TF_CT = tf(K_mixsyn_TF_b, K_mixsyn_TF_a);
K_mixsyn_TF_DT = c2d(K_mixsyn_TF_CT, Ts, 'zoh');

%% ---------------------------Data-driven----------------------------------%%
w = logspace(-2,2.55,50);
W2_frd = frd(W2,w);

% Initial stabilising controller
K_init = tf(0.001, 1, Ts);

% Frequency Grid
omegas = unique([logspace(log10 (0.3) , log10 (W2_frd.Frequency(end)) , 200) , ...
                 linspace(1, 10, 101), ...
                 linspace(100, 380, 101)]);

% LFR Models
G = augw(G_nominal, W1, W3, W2);
G = mktito(G, 1, ninputs);

% Controller Object 
K = datadriven.Controller.SS(ncontrolled_states, 1, ninputs, Ts);
K.setinitial(K_init);

% Synthesiser Object
synthesiser = Synthesiser(K);

% Add Hinf objectives
synthesiser.add_Hinf_objective(G, omegas);

% Add H2 objectives
% synthesiser.add_H2_objective(G, omegas);

% Add Hinf constraints
% synthesiser.add_Hinf_constraint(G, omegas);

% Add H2 objectives
% synthesiser.add_H2_constraint(G, omegas);

% Ensure that controller is stable
synthesiser.ensure_controller_stability(omegas);

% Optimize
output = synthesiser.synthesise();
K_DD = output.Controller;

%%
S = feedback(1, G_nominal * K_DD);
T = feedback(G_nominal * K_DD, 1);
U = feedback(K_DD, G_nominal);

figure;
step(T);
title('Step Response - Output (T) with W3');
grid on;

figure;
step(U);
title('Step Response - Control Signal (U) with W3');
grid on;

figure;
bodemag(S,1/W1);
legend('Sensitivity S(z)','$W1(z)^{-1}$','Interpreter','latex', 'Location', 'southeast');
title('Bode magnitude of sensitivity TF with W3');
grid on;

figure;
bodemag(T,1/W2);
legend('T(z)','$W2(z)^{-1}$','Interpreter','latex', 'Location', 'southeast');
title('Bode magnitude of T with W3');
grid on;

figure;
bodemag(U,1/W3);
legend('Input sensitivity U(z)','$W3(z)^{-1}$','Interpreter','latex', 'Location', 'southeast');
title('Bode magnitude of input sensitivity with W3');
grid on;


[K_DD_TF_b,K_DD_TF_a] = ss2tf(K_DD.A, K_DD.B, K_DD.C, K_DD.D);
K_DD_TF_CT = tf(K_DD_TF_b, K_DD_TF_a);
K_DD_TF_DT = c2d(K_DD_TF_CT, Ts, 'zoh');

