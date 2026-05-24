function plots_EDMD(n, Ts, mode, comparison, data_EDMD)
    % PLOTS
    
    % 0. Check folder
    if ~isfolder("figures")
        mkdir("figures");
    end
    
    % 1. States evolution (Trajectory 1) for the different models
    t = 0:Ts:(n.steps-1)*Ts;
    figure;
    hold on
    plot(t, comparison.q_nonlinear(1,:,1));
    plot(t, comparison.q_BILINEAR_DT(2,:,1)); %!!! Attention 2 psk 1 est la constante pour u
    plot(t, comparison.q_BILINEAR_CT(2,:,1));
    plot(t, comparison.q_EDMD(2,:,1));
    plot(t, comparison.q_LINEAR(1,:,1));
    title("State q(1) : rotary arm angle [rad]");
    legend("Nonlinear", "Bilinear DT", "Bilinear CT", "EDMD", "Linear");
    savefig("figures\q1.fig");
    saveas(gcf, 'figures\q1.png');
    hold off
    
    figure;
    hold on
    plot(t, comparison.q_nonlinear(2,:,1));
    plot(t, comparison.q_BILINEAR_DT(3,:,1));
    plot(t, comparison.q_BILINEAR_CT(3,:,1));
    plot(t, comparison.q_EDMD(3,:,1));
    plot(t, comparison.q_LINEAR(2,:,1));
    title("State q(2) : pendulum angle [rad]");
    legend("Nonlinear", "Bilinear DT", "Bilinear CT", "EDMD", "Linear");
    savefig("figures\q2.fig");
    saveas(gcf, 'figures\q2.png');
    hold off
    
    figure;
    hold on
    plot(t, comparison.q_nonlinear(3,:,1));
    plot(t, comparison.q_BILINEAR_DT(4,:,1));
    plot(t, comparison.q_BILINEAR_CT(4,:,1));
    plot(t, comparison.q_EDMD(4,:,1));
    plot(t, comparison.q_LINEAR(3,:,1));
    title("State q(3) : rotary arm speed [rad/s]");
    legend("Nonlinear", "Bilinear DT", "Bilinear CT", "EDMD", "Linear");
    savefig("figures\q3.fig");
    saveas(gcf, 'figures\q3.png');
    hold off
    
    figure;
    hold on
    plot(t, comparison.q_nonlinear(4,:,1));
    plot(t, comparison.q_BILINEAR_DT(5,:,1));
    plot(t, comparison.q_BILINEAR_CT(5,:,1));
    plot(t, comparison.q_EDMD(5,:,1));
    plot(t, comparison.q_LINEAR(4,:,1));
    title("State q(4) : pendulum speed [rad/s]");
    legend("Nonlinear", "Bilinear DT", "Bilinear CT", "EDMD", "Linear");
    savefig("figures\q4.fig");
    saveas(gcf, 'figures\q4.png');
    hold off
    
    % 2. Histograms of initial q2 state distribution
    figure;
    histogram(data_EDMD.training.q(2,1,:),BinWidth=2*pi/n.regions,BinLimits=[-pi,pi]);
    title("Training trajectories : distribution of q(2) at t=0 ");
    savefig("figures\q2_distrib_training_trajs.fig");
    saveas(gcf, 'figures\q2_distrib_training_trajs.png');
    
    
    figure;
    histogram(data_EDMD.testing.q(2,1,:),BinWidth=2*pi/n.regions,BinLimits=[-pi,pi]);
    title("Testing trajectories : distribution of q(2) at t=0 ");
    savefig("figures\q2_distrib_testing_trajs.fig");
    saveas(gcf, 'figures\q2_distrib_testing_trajs.png');
    
    
    % All q2 values across all time steps and trajectories
    q2_all_training = reshape(data_EDMD.training.q(2,:,:), 1, []);
    q2_all_testing  = reshape(data_EDMD.testing.q(2,:,:),  1, []);
    
    figure;
    % histogram(q2_all_training, BinWidth=2*pi/n.regions, BinLimits=[-pi,pi]);
    histfit(squeeze(q2_all_training),16);
    x = squeeze(q2_all_training);
    x = x(:);
    pd = fitdist(x, 'Normal');
    title('Training data : distribution of q_2 across all time steps');
    legend(sprintf('Variance : %d', pd.sigma));
    xlabel('\theta_2 (rad)'); ylabel('count');
    savefig('figures\q2_alltime_training.fig');
    saveas(gcf, 'figures\q2_alltime_training.png');
    
    
    figure;
    % histogram(q2_all_testing, BinWidth=2*pi/n.regions, BinLimits=[-pi,pi]);
    histfit(squeeze(q2_all_testing),16);
    x = squeeze(q2_all_testing);
    x = x(:);
    pd = fitdist(x, 'Normal');
    legend(sprintf('Variance : %d', pd.sigma));
    title('Testing data : distribution of q_2 across all time steps');
    xlabel('\theta_2 (rad)'); ylabel('count');
    savefig('figures\q2_alltime_testing.fig');
    saveas(gcf, 'figures\q2_alltime_testing.png');
    
    
    % 3. Bar charts of errors
    figure;
    x = ["LINEAR" "EDMD" "BILINEAR"];
    y = [comparison.errors.LINEAR_fro_training comparison.errors.EDMD_fro_training comparison.errors.BILINEAR_DT_fro_training];
    bar(x,y);
    title("Training errors (Frobenius norm)");
    savefig("figures\training_error_fro.fig");
    saveas(gcf, 'figures\training_error_fro.png');
    
    figure;
    x = ["EDMD" "BILINEAR"];
    y = [comparison.errors.EDMD_fro_training_lifted comparison.errors.BILINEAR_DT_fro_training_lifted];
    bar(x,y);
    title("Training errors in lifted space (Frobenius norm)");
    savefig("figures\training_error_fro_lifted.fig");
    saveas(gcf, 'figures\training_error_fro_lifted.png');
    
    
    figure;
    x = ["LINEAR" "EDMD" "BILINEAR"];
    y = [comparison.errors.LINEAR_fro_testing comparison.errors.EDMD_fro_testing comparison.errors.BILINEAR_DT_fro_testing];
    bar(x,y);
    title("Testing errors for original states (fro)");
    savefig("figures\testing_error_original_states_fro.fig");
    saveas(gcf, 'figures\testing_error_original_states_fro.png');
    
    
    figure;
    x = ["EDMD" "BILINEAR"];
    y = [comparison.errors.EDMD_fro_lifted_testing comparison.errors.BILINEAR_DT_fro_lifted_testing];
    bar(x,y);
    title("Testing errors for lifted states (fro)");
    savefig("figures\testing_error_lifted_states_fro.fig");
    saveas(gcf, 'figures\testing_error_lifted_states_fro.png');
    
    
    % 4. Phase portrait of q2 vs q4 (pendulum angle vs pendulum speed)
    q2_reshaped = reshape(data_EDMD.training.q(2,:,:), 1, []);
    q4_reshaped = reshape(data_EDMD.training.q(4,:,:), 1, []);
    figure;
    plot(q2_reshaped, q4_reshaped, '--')
    xlabel('q2 [rad]'); ylabel('q4 [rad/s]');
    title('Phase portrait of q2 vs q4');
    savefig("figures\phase_portrait_q2_vs_q4.fig");
    saveas(gcf, 'figures\phase_portrait_q2_vs_q4.png');
    
    
    % 5. Phase portrait of q1 vs q3 (rotary angle vs rotary speed)
    q1_reshaped = reshape(data_EDMD.training.q(1,:,:), 1, []);
    q3_reshaped = reshape(data_EDMD.training.q(3,:,:), 1, []);
    figure;
    plot(q1_reshaped, q3_reshaped, '--')
    xlabel('q1 [rad]'); ylabel('q3 [rad/s]');
    title('Phase portrait of q1 vs q3');
    savefig("figures\phase_portrait_q1_vs_q3.fig");
    saveas(gcf, 'figures\phase_portrait_q1_vs_q3.png');
    
    
    
    % 6. Input trajectories (to see if suficiently excited)
    figure;
    hold on;
    for i = 1:n.trajs_training
        plot(data_EDMD.training.u(:,:,i));
    end
    yline(15, 'r--', 'n.umax');
    yline(-15, 'r--', '-n.umax');
    yline(n.umax, 'r--', 'nmax');
    yline(-n.umax, 'r--', '-nmax');
    title('Training trajectories : input u');
    xlabel('Steps');
    ylabel('u [V]');
    ylim([-16 16]);
    hold off;
    savefig("figures\u_training.fig");
    saveas(gcf, 'figures\u_training.png');
    
    % 7. Singular Values of Phi for the models (to see if some lifting functions are poorly chosen)
    figure;
    semilogy(comparison.svds.LINEAR);
    title('Singular values of lifted states matrix for linear model');
    savefig("figures\svd_linear.fig");
    saveas(gcf, 'figures\svd_linear.png');
    
    
    figure;
    semilogy(comparison.svds.EDMD);
    title('Singular values of lifted states matrix for EDMD linear model');
    savefig("figures\svd_edmd.fig");
    saveas(gcf, 'figures\svd_edmd.png');
    
    
    figure;
    semilogy(comparison.svds.BILINEAR_CT);
    title('Singular values of lifted states matrix for EDMD bilinear model');
    savefig("figures\svd_bilinear.fig");
    saveas(gcf, 'figures\svd_bilinear.png');
    
    %%% CONFIG FILE TO KNOW ALL PARAMETERS USED
    
    timestamp = datetime('now','Format','yyyy_MM_dd_HHmmss');
    filename = "figures\report_" + string(timestamp) + ".txt";
    f = fopen(filename,'w');
    
    fprintf(f,"==== TEST CONFIGURATION ====\n");
    fprintf(f,"Ts = %.3f\n", Ts);
    fprintf(f,"mode = %s\n", mode);
    fprintf(f,"n.umax = %.3f\n", n.umax);
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
    fprintf(f,"BILINEAR Fro = %.3f\n", comparison.errors.BILINEAR_DT_fro_training);
    % fprintf(f,"LINEAR 2norm = %.3f\n", comparison.errors.LINEAR_2norm_training);
    % fprintf(f,"EDMD 2norm = %.3f\n", comparison.errors.EDMD_2norm_training);
    % fprintf(f,"BILINEAR 2norm = %.3f\n", comparison.errors.BILINEAR_2norm_training);
    
    fprintf(f,"\n--- TESTING ERRORS (true states) ---\n");
    fprintf(f,"LINEAR fro = %.3f\n", comparison.errors.LINEAR_fro_testing);
    fprintf(f,"EDMD fro = %.3f\n", comparison.errors.EDMD_fro_testing);
    fprintf(f,"BILINEAR fro = %.3f\n", comparison.errors.BILINEAR_DT_fro_testing);
    
    fprintf(f,"\n--- TESTING ERRORS (lifted states) ---\n");
    fprintf(f,"EDMD lifted fro = %.3f\n", comparison.errors.EDMD_fro_lifted_testing);
    fprintf(f,"BILINEAR lifted fro = %.3f\n", comparison.errors.BILINEAR_DT_fro_lifted_testing);
    
    fprintf(f,"\n--- Condition number of lifted_Q ---\n");
    fprintf(f,"LINEAR : %.3f\n", comparison.cond.LINEAR);
    fprintf(f,"EDMD : %.3f\n", comparison.cond.EDMD);
    fprintf(f,"BILINEAR : %.3f\n", comparison.cond.BILINEAR_CT);
    
    fclose(f);
    
    disp("Report saved: " + filename);



end