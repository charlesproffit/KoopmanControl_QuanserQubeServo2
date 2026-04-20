clc;
clear;
close all;
% addpath("D:\MA4\PDSE_II\MatLab\LiftingFunctions");

rng(22);

% System Parameters
Ts = 0.008;         % Time step
n.states = 4;       % Number of states
n.inputs = 1;       % Number of inputs   
n.outputs = 1;      % Number of outputs
n.regions = 16;     % Number of subdivisions of a full circle for the initial conditions

n.train_test_ratio = 0.8;
n.cross_val_groups = 5;
umax = 6;          % Max possible input
mode = "REAL";

if mode == "SIM"
    % Experiment Parameters
    n.trajs = 100;          % Number of trajectories generated
    n.steps = 500;          % Total number of steps per trajectory (u(k)=0 if k>n.steps_excited)
end

[f_continuous, f_discrete] = dynamics(Ts);
[SS_LINEAR_CT, SS_LINEAR_DT] = dynamics_linearized([0;0;0;0],0,Ts);
[K_LQR,~,~] = dlqr(SS_LINEAR_DT.A, SS_LINEAR_DT.B, diag([0.41, 0, 0, 0]), 0.016);

% Koopman functions
f_lifting = @(q)[
    q(1);
    q(2);
    q(3);
    q(4);
    sin(q(2));
    cos(q(2));
    sin(q(2))^2;
    cos(q(2))*sin(q(2))*q(3)*q(4);
    sin(q(2))*q(4)^2;
    cos(q(2))*sin(q(2))*q(3)^2;
    % sin(q(2))*q(3);
    % cos(q(2))*q(3);
    % cos(q(2))*q(4)
];
f_linear = @(q)[q(1); q(2); q(3); q(4)];
% f_lifting = f_linear;
n.lifted_states = size(f_lifting(zeros(n.states,1)),1);

%--------------DATA COLLECTION-----------------------%

% Data collected with non-linear ODEs model
if mode == "REAL"
    [data_EDMD,n] = structure_data('data.mat', n, 'all');
else
    [data_EDMD,n] = collect_data(f_discrete, umax, n, 'EDMD', K_LQR);
end

%--------------STATE-SPACE IDENTIFICATION WITH EDMD-----------------------%

% Bilinear EDMD lifted model
[M_BILINEAR, ~] = compute_edmd(data_EDMD, f_lifting, n, 'BILINEAR', 'RIDGE');

% Linear EDMD lifted model
[M_EDMD, ~] = compute_edmd(data_EDMD, f_lifting, n, 'EDMD', 'RIDGE');

% Linear EDMD non-lifted model
[M_LINEAR, ~] = compute_edmd(data_EDMD, f_linear, n, 'LINEAR', 'RIDGE');

% Comparison of the 4 models
comparison = compare_models(data_EDMD, f_lifting, M_BILINEAR, M_EDMD, M_LINEAR, n);

% plots(n, Ts, mode, umax, comparison, data_EDMD, M_LINEAR, M_EDMD, M_BILINEAR);


%--------------DATA-DRIVEN FEEDBACK LINEARIZATION AND CONTROL-----------------------%
M_DDFL = feedback_linearization(f_lifting, M_BILINEAR, n, 1);

K = dlqr(M_DDFL.A, M_DDFL.B, eye(size(M_DDFL.A, 1)), 1);

% Reference
q_ref = [0;0;0;0];
z_ref = M_DDFL.T(q_ref);

trajs = simul_control(f_discrete, umax, n, K_LQR, M_DDFL, K, z_ref);

%---------FREQUENCY RESPONSE IDENTIFICATION WITH FOURIER ANALYSIS---------%

% % Data collected with PRBS signal
% data_FRF = collect_data(f_discrete, umax, n, 'FRF');
% 
% % FRF model
% [M_FRF, data_FRF_lifted] = compute_frf(data_FRF, f_lifting, Ts, n, 'FRF');
% 
% 
% opts = bodeoptions('cstprefs');
% opts.PhaseWrapping='on';
% bode(M_FRF, '.', ss(SS_LINEAR_CT.A, SS_LINEAR_CT.B, [0,1,0,0], []), 'r.', opts)
% compare(M_FRF,ss(SS_LINEAR_CT.A, SS_LINEAR_CT.B, [0,1,0,0], []))
