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
f_discrete_with_LQR = @(q,u) f_discrete(q,u - K_LQR*q);

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
];
f_linear = @(q)[q(1); q(2); q(3); q(4)];
n.lifted_states = size(f_lifting(zeros(n.states,1)),1);

%--------------DATA COLLECTION-----------------------%

% Data collected with non-linear ODEs model
if mode == "REAL"
    [data_EDMD,n] = structure_data('data.mat', n, 'all');
else
    [data_EDMD,n] = collect_data(f_discrete_with_LQR, umax, n, 'EDMD');
end

%--------------STATE-SPACE IDENTIFICATION WITH EDMD-----------------------%

% Bilinear EDMD lifted model
[M_BILINEAR_CT, M_BILINEAR_DT, ~] = compute_edmd(data_EDMD, f_lifting, n, Ts, 'BILINEAR', 'RIDGE');

% Linear EDMD lifted model
[~, M_EDMD_DT, ~] = compute_edmd(data_EDMD, f_lifting, n, Ts, 'EDMD', 'RIDGE');

% Linear EDMD non-lifted model
[~, M_LINEAR_DT, ~] = compute_edmd(data_EDMD, f_linear, n, Ts, 'LINEAR', 'RIDGE');

% Comparison of the 4 models
comparison = compare_models(data_EDMD, f_lifting, M_BILINEAR_DT, M_EDMD_DT, M_LINEAR_DT, n);

% plots_EDMD(n, Ts, mode, umax, comparison, data_EDMD, M_LINEAR_DT, M_EDMD_DT, M_BILINEAR_DT);


%--------------DATA-DRIVEN FEEDBACK LINEARIZATION AND CONTROL-----------------------%
[M_DDFL_CT, M_DDFL_DT] = feedback_linearization(f_lifting, M_BILINEAR_CT, n, 1, Ts);

K = lqr(M_DDFL_CT.A, M_DDFL_CT.B, 0.01*eye(size(M_DDFL_CT.A, 1)), 0.016);

% Reference
q_ref = [0;0;0;0];
z_ref = M_DDFL_DT.T(q_ref);

trajs = simul_control(f_discrete_with_LQR, umax, n, M_DDFL_DT, K, z_ref);