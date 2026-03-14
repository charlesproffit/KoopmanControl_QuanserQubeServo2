clc;
clear;
close all;
addpath("D:\MA4\PDSE_II\MatLab\LiftingFunctions");

% System Parameters
Ts = 0.01;          % Time step
f_linear = @(q)[q(1); q(2); q(3); q(4)];
n.states = 4;       % Number of states
n.inputs = 1;       % Number of inputs   
n.outputs = 1;      % Number of outputs
umax = 10;          % Max possible input
mode = "SIM";

[f_continuous, f_discrete] = dynamics(Ts);

% Test linear f_discrete to see if EDMD works
A = [0.5 0 0.1 0;
     0 0.6 0 0.1;
     0.2 0 0.7 0;
     0 0.1 0 0.8];
B = [0;
     1;
     0;
     0];
f_discrete = @(q,u)A*f_linear(q)+B*u;

% Koopman Functions Parameters
f_lifting = @(q)[
    q(1);
    q(2);
    q(3);
    q(4);
    cos(q(2))*sin(q(2))*q(3)*q(4);
    % sin(q(2))*q(4)^2;
    % cos(q(2))*sin(q(2))*q(3)^2;
    % sin(q(2))*cos(q(2));
    % sin(q(2));
    % q(3)^2;
    % q(4)^2;
    % q(3)*q(4);
    % sin(q(2))*q(3);
    % sin(q(2))*q(4);
    % cos(q(2))*q(3);
    % cos(q(2))*q(4)
];
f_lifting = f_linear;
n.lifted_states = size(f_lifting(zeros(n.states,1)),1);

%--------------STATE-SPACE IDENTIFICATION WITH EDMD-----------------------%

% Data collected with non-linear ODEs model
if mode == "REAL"
    n.train_test_ratio = 0.6;
    n.regions = 16;         % Number of subdivisions of a full circle for the initial conditions
    [data_EDMD,n] = structure_data("data.tdms", n);
else
    % Experiment Parameters
    n.regions = 16;         % Number of subdivisions of a full circle for the initial conditions
    n.trajs = 100;          % Number of trajectories generated
    n.steps = 100;          % Total number of steps per trajectory (u(k)=0 if k>n.steps_excited)
    n.train_test_ratio = 0.6;
    [data_EDMD,n] = collect_data(f_discrete, umax, n, 'EDMD');
end

% Bilinear EDMD lifted model
[M_BILINEAR, ~] = compute_edmd(data_EDMD, f_lifting, n, 'BILINEAR');

% Linear EDMD lifted model
[M_EDMD, ~] = compute_edmd(data_EDMD, f_lifting, n, 'EDMD');

% Linear EDMD non-lifted model
[M_LINEAR, ~] = compute_edmd(data_EDMD, f_linear, n, 'LINEAR');

% Comparison of the 4 models
comparison = compare_models(data_EDMD, f_lifting, M_BILINEAR, M_EDMD, M_LINEAR, n);

t = 0:Ts:(n.steps)*Ts;
figure;
hold on
plot(t, comparison.q_nonlinear(1,:,1));
plot(t, comparison.q_BILINEAR(1,:,1));
plot(t, comparison.q_EDMD(1,:,1));
plot(t, comparison.q_LINEAR(1,:,1));
title("State q(1)");
legend("Nonlinear", "Bilinear", "EDMD", "Linear");
hold off

figure;
hold on
plot(t, comparison.q_nonlinear(2,:,1));
plot(t, comparison.q_BILINEAR(2,:,1));
plot(t, comparison.q_EDMD(2,:,1));
plot(t, comparison.q_LINEAR(2,:,1));
title("State q(2)");
legend("Nonlinear", "Bilinear", "EDMD", "Linear");
hold off

figure;
hold on
plot(t, comparison.q_nonlinear(3,:,1));
plot(t, comparison.q_BILINEAR(3,:,1));
plot(t, comparison.q_EDMD(3,:,1));
plot(t, comparison.q_LINEAR(3,:,1));
title("State q(3)");
legend("Nonlinear", "Bilinear", "EDMD", "Linear");
hold off

figure;
hold on
plot(t, comparison.q_nonlinear(4,:,1));
plot(t, comparison.q_BILINEAR(4,:,1));
plot(t, comparison.q_EDMD(4,:,1));
plot(t, comparison.q_LINEAR(4,:,1));
title("State q(4)");
legend("Nonlinear", "Bilinear", "EDMD", "Linear");
hold off

%---------FREQUENCY RESPONSE IDENTIFICATION WITH FOURIER ANALYSIS---------%

% % Data collected with PRBS signal
% data_FRF = collect_data(f_discrete, umax, n, 'FRF');
% 
% % FRF model
% [M_FRF, data_FRF_lifted] = compute_frf(data_FRF, f_lifting, Ts, n, 'FRF');
% 
% [SS_LINEAR_CT, SS_LINEAR_DT] = dynamics_linearized([0;0;0;0],0,Ts);
% 
% opts = bodeoptions('cstprefs');
% opts.PhaseWrapping='on';
% bode(M_FRF, '.', ss(SS_LINEAR_CT.A, SS_LINEAR_CT.B, [0,1,0,0], []), 'r.', opts)
% compare(M_FRF,ss(SS_LINEAR_CT.A, SS_LINEAR_CT.B, [0,1,0,0], []))