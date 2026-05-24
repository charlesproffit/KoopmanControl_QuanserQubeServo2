clc;
clear;
close all;
addpath(genpath('lib'));
addpath(genpath('data'));

rng(43);

    % States q :
    %   q(1) => Angle of rotary arm
    %   q(2) => Angle of pendulum arm
    %   q(3) => Angular speed of rotary arm q(1)_dot
    %   q(4) => Angular speed of pendulum arm q(2)_dot

%% ------------- SYSTEM PARAMETERS ---------------------------- %%
Ts = 0.008;         % Time step
n.states = 4;       % Number of states
n.inputs = 1;       % Number of inputs   
n.outputs = 1;      % Number of outputs
n.regions = 16;     % Number of subdivisions of a full circle for the initial conditions
n.controlled_states = 2;
n.train_test_ratio = 0.8;
n.cross_val_groups = 5;
n.trajs_control = 2;
n.umax = 15;          % Max possible input
mode = "SIM";

if mode == "SIM"
    % Experiment Parameters
    n.trajs = 500;          % Number of trajectories generated
    n.steps = 500;          % Total number of steps per trajectory (u(k)=0 if k>n.steps_excited)
    n.trajs_training = fix(n.train_test_ratio*n.trajs);
    n.trajs_testing = n.trajs - n.trajs_training;
end

%% -------- DYNAMICS DEFINITION --------------------------- %%
% Design of internal LQR
[SS_LINEAR_CT, SS_LINEAR_DT] = dynamics_linearized([0;0;0;0],0,Ts);
[K_LQR, ~, ~] = dlqr(SS_LINEAR_DT.A, SS_LINEAR_DT.B, diag([0.41, 0, 0, 0]), 0.016);
% K_LQR = [4,1,0,0];
K_LQR = [0,0,0,0];

% Dynamics
[f_continuous, f_discrete, funcs] = dynamics(Ts, 'euler');
f_discrete_with_LQR = @(q,u) f_discrete(q,u - K_LQR*q);

% f_lifting = @(q)[
%     1;
%     q(1); q(2); q(3); q(4);
%     sin(q(2));
%     cos(q(2));
%     sin(q(2))*cos(q(2));
%     sin(q(2))^2;
%     q(3)^2;
%     q(4)^2;
%     q(3)*q(4);
%     sin(q(2))*q(4)^2;
%     cos(q(2))*q(4);
%     sin(q(2))*cos(q(2))*q(3)*q(4);
%     sin(q(2))*cos(q(2))*q(3)^2;
%     sin(q(2))*cos(q(2))^2*q(3)^2;
%     % 1 / funcs.D(q);
%     % sin(q(2)) / funcs.D(q);
%     % cos(q(2)) / funcs.D(q);
%     % sin(q(2))^2 / funcs.D(q);
%     % q(3) / funcs.D(q);
% ];

f_lifting = @(q)[
    1;
    q;
    % funcs.f(q);
    funcs.f_not_redundant(q);
    % funcs.f(q) - funcs.g(q)*K_LQR*q;
    funcs.g_not_redundant(q);
    % funcs.g(q);
];
f_linear = @(q)[q(1); q(2); q(3); q(4)];
n.lifted_states = size(f_lifting(zeros(n.states,1)),1);

%% --------------DATA COLLECTION----------------------- %%

if mode == "REAL"
    [data_EDMD,n] = structure_data('data.mat', n, 'all');
else
    % Data collected with non-linear ODEs model
    [data_EDMD,n] = collect_data(f_discrete_with_LQR, n);
end

%% --------------STATE-SPACE IDENTIFICATION WITH EDMD----------------- %%

% Bilinear EDMD lifted model
[M_BILINEAR_CT, M_BILINEAR_DT] = compute_EDMD(data_EDMD, f_lifting, n, Ts, 'BILINEAR', 'LS');
% M_BILINEAR_CT.A = M_BILINEAR_CT.A + M_BILINEAR_CT.N*K_LQR*M_BILINEAR_CT.C;
save("data\M_BILINEAR_CT.mat", "M_BILINEAR_CT");

tol = 1e-9;  % adjust to taste
A_clean = M_BILINEAR_CT.A .* (abs(M_BILINEAR_CT.A) > tol);
N_clean = M_BILINEAR_CT.N .* (abs(M_BILINEAR_CT.N) > tol);

disp('A_clean:'); disp(A_clean)
disp('N_clean:'); disp(N_clean)

% Linear EDMD lifted model
[~, M_EDMD_DT] = compute_EDMD(data_EDMD, f_lifting, n, Ts, 'EDMD', 'LS');

% Linear EDMD non-lifted model
[~, M_LINEAR_DT] = compute_EDMD(data_EDMD, f_linear, n, Ts, 'LINEAR', 'LS');

% Comparison of the 4 models
comparison = compare_models_EDMD(data_EDMD, f_lifting, M_BILINEAR_DT, M_BILINEAR_CT, M_EDMD_DT, M_LINEAR_DT, n, Ts);

plots_EDMD(n, Ts, mode, comparison, data_EDMD);


%% ---------DATA-DRIVEN FEEDBACK LINEARIZATION SYSTEM-------------------- %%
M_DDFL_CT = compute_DDFL(f_lifting, M_BILINEAR_CT, n, 1e-6);
M_MBFL_CT = compute_MBFL(K_LQR);

% scaling
% x0 = [0; 0; 0; 0];
% gamma_D = M_DDFL_CT.gamma(x0);
% gamma_M = M_MBFL_CT.gamma(x0);
% scale = gamma_M / gamma_D;
% gamma_old = M_DDFL_CT.gamma;
% % etta_old  = M_DDFL_CT.etta;
% M_DDFL_CT.gamma = @(q)(scale * gamma_old(q));
% % M_DDFL_CT.etta  = @(q)(scale * etta_old(q));

save("data\M_MBFL_CT.mat", "M_MBFL_CT");
save("data\M_DDFL_CT.mat", "M_DDFL_CT");

% trajs_OL = compare_models_FL_OL(M_MBFL_CT, M_DDFL_CT, f_discrete_with_LQR, n, Ts);

%% ----------DATA-DRIVEN FEEDBACK LINEARIZATION CONTROLLER -------------------%%
Q = [
    1, 0;
    0, 1;
];
R = 0.01;
K = lqr(M_DDFL_CT.A, M_DDFL_CT.B, Q, R); % same A and B for MBFL and DDFL so valid for both simulations

q_ref = [0;0;0;0];

disp('Collecting trajectories CL, DDFL');
trajs_CL_DDFL = simul_control(f_discrete_with_LQR, n, M_DDFL_CT, K, q_ref);
disp('Collecting trajectories CL, MBFL');
trajs_CL_MBFL = simul_control(f_discrete_with_LQR, n, M_MBFL_CT, K, q_ref);
plots_FL_CL(trajs_CL_DDFL, Ts, 1000, 'CL_{DDFL}, SIM data');
plots_FL_CL(trajs_CL_MBFL, Ts, 1000, 'CL_{MBFL}, SIM data');
compare_models_FL_CL(trajs_CL_MBFL, trajs_CL_DDFL, 'CL_{MBFL} Vs CL_{DDFL}, SIM DATA');

%% Comparison with no controller
M_OL.etta = @(q)(0);
M_OL.gamma = @(q)(1);
M_OL.T = @(q)(q);
[K_OL, ~, ~] = dlqr(SS_LINEAR_DT.A, SS_LINEAR_DT.B, diag([0.41, 0, 0.01, 0]), 0.01);
% n.train_test_ratio = 1;
% n.umax = 0;
% [data_EDMD_2,n] = collect_data(f_discrete_with_LQR, n);
% trajs_OL_LQR_only = data_EDMD_2.training;
disp('Collecting trajectories with external K on linearized system');
trajs_OL_LQR_only = simul_control(f_discrete_with_LQR, n, M_OL, K_OL, q_ref);
plots_FL_CL(trajs_OL_LQR_only, Ts, 1000, 'OL, SIM data');
compare_models_FL_CL(trajs_OL_LQR_only, trajs_CL_MBFL, 'OL Vs CL_{MBFL}, SIM data'); %%Changed
compare_models_FL_CL(trajs_OL_LQR_only, trajs_CL_DDFL, 'OL Vs CL_{DDFL}, SIM data'); %%Changed


%% Comparison with real data DDFL model
[data_EDMD_real,n] = structure_data('data.mat', n, 'all');
[M_BILINEAR_CT_real, ~] = compute_EDMD(data_EDMD_real, f_lifting, n, Ts, 'BILINEAR', 'LS');
M_DDFL_CT_real = compute_DDFL(f_lifting, M_BILINEAR_CT_real, n, 1e-6);

% %scaling
% x0 = [0; 0; 0; 0];
% gamma_D = M_DDFL_CT_real.gamma(x0);
% gamma_M = M_MBFL_CT.gamma(x0);
% scale = gamma_M / gamma_D;
% gamma_old = M_DDFL_CT_real.gamma;
% M_DDFL_CT_real.gamma = @(q)(scale * gamma_old(q));

save("data\M_DDFL_CT_real.mat", "M_DDFL_CT_real");

n.steps = 500;
disp('Collecting trajectories CL, DDFL with real data');
trajs_CL_DDFL_real = simul_control(f_discrete_with_LQR, n, M_DDFL_CT_real, K, q_ref);
plots_FL_CL(trajs_CL_DDFL_real, Ts, 1000, 'CL_{DDFL}, REAL data');
compare_models_FL_CL(trajs_CL_DDFL, trajs_CL_DDFL_real, 'CL_{DDFL}, SIM data Vs CL_{DDFL}, REAL data');




%% -------------- ROBUST CONTROL ------------------------------ %%
% 
% 
% 





