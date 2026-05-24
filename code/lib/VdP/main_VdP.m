clc;
clear;
close all;

addpath(genpath('D:/MA4/PDSE_II/MatLab/main/lib'));
addpath(genpath('D:/MA4/PDSE_II/MatLab/main/data'));

rng(43);

    % States q :
    %   q(1) => Angle of rotary arm
    %   q(2) => Angle of pendulum arm
    %   q(3) => Angular speed of rotary arm q(1)_dot
    %   q(4) => Angular speed of pendulum arm q(2)_dot

%% --------------- TEST PIPELINE WITH KNOWN VDP MODEL ---------------- %%

Ts = 1e-2;         % Time step
n.states = 2;       % Number of states
n.inputs = 1;       % Number of inputs   
n.outputs = 1;      % Number of outputs
n.controlled_states = 2;

n.train_test_ratio = 0.8;
n.cross_val_groups = 5;
n.umax = 5;          % Max possible input
mode = "SIM";
n.trajs_control = 5;


if mode == "SIM"
    % Experiment Parameters
    n.trajs = 2;          % Number of trajectories generated
    n.steps = 2000;          % Total number of steps per trajectory (u(k)=0 if k>n.steps_excited)
    n.trajs_training = fix(n.train_test_ratio*n.trajs);
    n.trajs_testing = n.trajs - n.trajs_training;
end


[f_continuous_VdP, f_discrete_VdP] = dynamics_VdP(Ts, 'euler');

f_lifting = @(q)[ % Hermite Polynomials up to order 3
    1; % order 0
    2*q(1); 2*q(2); % order 1
    4*q(1)^2-2; 4*q(1)*q(2); 4*q(2)^2-2; % order 2
    8*q(1)^3-12*q(1); 8*q(1)^2*q(2)-12*q(2); 8*q(1)*q(2)^2-12*q(1); 8*q(2)^3-12*q(2); % order 3 - comment out for results in Fig. 3
];
f_linear = @(q)[q(1); q(2)];
n.lifted_states = size(f_lifting(zeros(n.states,1)),1);

[data_EDMD,n] = collect_data_VdP(f_discrete_VdP, n, 'EDMD');

% Bilinear EDMD lifted model
[M_BILINEAR_CT, M_BILINEAR_DT] = compute_EDMD(data_EDMD, f_lifting, n, Ts, 'BILINEAR', 'LS');
save("M_BILINEAR_CT_VdP.mat", "M_BILINEAR_CT");

% Linear EDMD lifted model
[~, M_EDMD_DT] = compute_EDMD(data_EDMD, f_lifting, n, Ts, 'EDMD', 'LS');

% Linear EDMD non-lifted model
[~, M_LINEAR_DT] = compute_EDMD(data_EDMD, f_linear, n, Ts, 'LINEAR', 'LS');

% Comparison of the 4 models
comparison = compare_models_EDMD(data_EDMD, f_lifting, M_BILINEAR_DT, M_BILINEAR_CT, M_EDMD_DT, M_LINEAR_DT, n, Ts);

% plots_EDMD(n, Ts, mode, comparison, data_EDMD, M_LINEAR_DT, M_EDMD_DT, M_BILINEAR_DT, M_BILINEAR_CT);

M_DDFL_CT = compute_DDFL_VdP(f_lifting, M_BILINEAR_CT, n, 1);
M_MBFL_CT = compute_MBFL_VdP();
trajs = compare_models_FL_VdP(M_MBFL_CT, M_DDFL_CT, f_discrete_VdP, n, Ts);
plots_FL_OL(trajs, Ts, 1000);

K = [0,0];
q_ref = [0;0];
trajs_CL_DDFL = simul_control_VdP(f_discrete_VdP, n, M_DDFL_CT, K, q_ref);
plots_FL_CL(trajs_CL_DDFL, Ts, 1000, 'CL_{DDFL}, REAL data');
trajs_CL_MBFL = simul_control_VdP(f_discrete_VdP, n, M_MBFL_CT, K, q_ref);
plots_FL_CL(trajs_CL_MBFL, Ts, 1000, 'CL_{MBFL}, REAL data');



