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
n.regions = 16;     % Number of subdivisions of a full circle for the initial conditions

n.train_test_ratio = 0.6;
umax = 10;          % Max possible input
mode = "REAL";

[f_continuous, f_discrete] = dynamics(Ts);

% Test linear f_discrete to see if EDMD works
% A = [0.5 0 0.1 0;
%      0 0.6 0 0.1;
%      0.2 0 0.7 0;
%      0 0.1 0 0.8];
% B = [0;
%      1;
%      0;
%      0];
% f_discrete = @(q,u)A*f_linear(q)+B*u;

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
% f_lifting = f_linear;
n.lifted_states = size(f_lifting(zeros(n.states,1)),1);

%--------------STATE-SPACE IDENTIFICATION WITH EDMD-----------------------%

% Data collected with non-linear ODEs model
if mode == "REAL"
    [data_EDMD,n] = structure_data("data.tdms", n);
else
    % Experiment Parameters
    n.trajs = 100;          % Number of trajectories generated
    n.steps = 100;          % Total number of steps per trajectory (u(k)=0 if k>n.steps_excited)
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

% PLOTS

% 1. States evolution (Trajectory 1) for the different models
t = 0:Ts:(n.steps)*Ts;
figure;
hold on
plot(t, comparison.q_nonlinear(1,:,1));
plot(t, comparison.q_BILINEAR(1,:,1));
plot(t, comparison.q_EDMD(1,:,1));
plot(t, comparison.q_LINEAR(1,:,1));
title("State q(1)");
legend("Nonlinear", "Bilinear", "EDMD", "Linear");
savefig("q1.fig");
hold off

figure;
hold on
plot(t, comparison.q_nonlinear(2,:,1));
plot(t, comparison.q_BILINEAR(2,:,1));
plot(t, comparison.q_EDMD(2,:,1));
plot(t, comparison.q_LINEAR(2,:,1));
title("State q(2)");
legend("Nonlinear", "Bilinear", "EDMD", "Linear");
savefig("q2.fig");
hold off

figure;
hold on
plot(t, comparison.q_nonlinear(3,:,1));
plot(t, comparison.q_BILINEAR(3,:,1));
plot(t, comparison.q_EDMD(3,:,1));
plot(t, comparison.q_LINEAR(3,:,1));
title("State q(3)");
legend("Nonlinear", "Bilinear", "EDMD", "Linear");
savefig("q3.fig");
hold off

figure;
hold on
plot(t, comparison.q_nonlinear(4,:,1));
plot(t, comparison.q_BILINEAR(4,:,1));
plot(t, comparison.q_EDMD(4,:,1));
plot(t, comparison.q_LINEAR(4,:,1));
title("State q(4)");
legend("Nonlinear", "Bilinear", "EDMD", "Linear");
savefig("q4.fig");
hold off

% 2. Histograms of initial q2 state distribution
figure;
data_EDMD.training.q(2,1,:)
histogram(data_EDMD.training.q(2,1,:),BinWidth=2*pi/n.regions,BinLimits=[-pi,pi]);
title("Training trajectories : distribution of q(2) at t=0 ");
savefig("q2_distrib_training_trajs.fig");


figure;
histogram(data_EDMD.testing.q(2,1,:),BinWidth=2*pi/n.regions,BinLimits=[-pi,pi]);
title("Testing trajectories : distribution of q(2) at t=0 ");
savefig("q2_distrib_testing_trajs.fig");

% 3. Bar charts of errors
figure;
x = ["LINEAR" "EDMD" "BILINEAR"];
title("Training errors (Frobenius norm)");
y = [comparison.errors.LINEAR_fro_training comparison.errors.EDMD_fro_training comparison.errors.BILINEAR_fro_training];
bar(x,y);
savefig("training_error_fro.fig");


figure;
x = ["LINEAR" "EDMD" "BILINEAR"];
title("Training errors (2-norm)");
y = [comparison.errors.LINEAR_2norm_training comparison.errors.EDMD_2norm_training comparison.errors.BILINEAR_2norm_training];
bar(x,y);
savefig("training_error_2norm.fig");

figure;
x = ["LINEAR" "EDMD" "BILINEAR"];
title("Testing errors for original states (2-norm)");
y = [comparison.errors.LINEAR_2norm comparison.errors.EDMD_2norm comparison.errors.BILINEAR_2norm];
bar(x,y);
savefig("testing_error_original_states_2norm.fig");

figure;
x = ["EDMD" "BILINEAR"];
title("Testing errors for lifted states (2-norm)");
y = [comparison.errors.EDMD_2norm_lifted comparison.errors.BILINEAR_2norm_lifted];
bar(x,y);
savefig("testing_error_lifted_states_2norm.fig");


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


%%% CONFIG FILE TO KNOW ALL PARAMETERS USED
timestamp = datetime('now','Format','yyyy_MM_dd_HHmmss');
filename = "report_" + string(timestamp) + ".txt";
f = fopen(filename,'w');

fprintf(f,"==== TEST CONFIGURATION ====\n");
fprintf(f,"Ts = %.3f\n", Ts);
fprintf(f,"mode = %s\n", mode);
fprintf(f,"umax = %.3f\n", umax);
fprintf(f,"\n--- n parameters ---\n");
fprintf(f,"n.states = %d\n", n.states);
fprintf(f,"n.inputs = %d\n", n.inputs);
fprintf(f,"n.outputs = %d\n", n.outputs);
fprintf(f,"n.lifted_states = %d\n", n.lifted_states);
fprintf(f,"n.trajs = %d\n", n.trajs);
fprintf(f,"n.steps = %d\n", n.steps);
fprintf(f,"n.train_test_ratio = %.3f\n", n.train_test_ratio);
fprintf(f,"n.trajs_training = %d\n", n.trajs_training);
fprintf(f,"n.trajs_testing = %d\n", n.trajs_testing);
fprintf(f,"n.lifted_states = %d\n", n.lifted_states);

fprintf(f,"\n--- TRAINING ERRORS ---\n");
fprintf(f,"LINEAR Fro = %.3f\n", comparison.errors.LINEAR_fro_training);
fprintf(f,"EDMD Fro = %.3f\n", comparison.errors.EDMD_fro_training);
fprintf(f,"BILINEAR Fro = %.3f\n", comparison.errors.BILINEAR_fro_training);
fprintf(f,"LINEAR 2norm = %.3f\n", comparison.errors.LINEAR_2norm_training);
fprintf(f,"EDMD 2norm = %.3f\n", comparison.errors.EDMD_2norm_training);
fprintf(f,"BILINEAR 2norm = %.3f\n", comparison.errors.BILINEAR_2norm_training);

fprintf(f,"\n--- TESTING ERRORS (true states) ---\n");
fprintf(f,"LINEAR 2norm = %.3f\n", comparison.errors.LINEAR_2norm);
fprintf(f,"EDMD 2norm = %.3f\n", comparison.errors.EDMD_2norm);
fprintf(f,"BILINEAR 2norm = %.3f\n", comparison.errors.BILINEAR_2norm);

fprintf(f,"\n--- TESTING ERRORS (lifted states) ---\n");
fprintf(f,"EDMD lifted 2norm = %.3f\n", comparison.errors.EDMD_2norm_lifted);
fprintf(f,"BILINEAR lifted 2norm = %.3f\n", comparison.errors.BILINEAR_2norm_lifted);

fclose(f);

disp("Report saved: " + filename);