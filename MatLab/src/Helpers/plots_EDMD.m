function plots_EDMD(Ts, comparison)
    % PLOTS
    nsteps = size(comparison.X_nonlinear,2);
    ntrajs_training = size(comparison.X_training,3);
    
    % 0. Check folder
    if ~isfolder("figures")
        mkdir("figures");
    end
    
    % 1. States evolution (Trajectory 1) for the different models
    t = 0:Ts:(nsteps-1)*Ts;
    figure;
    hold on
    plot(t, comparison.X_nonlinear(1,:,1));
    plot(t, comparison.X_BILIN(1,:,1));
    plot(t, comparison.X_EDMD(1,:,1));
    plot(t, comparison.X_LIN(1,:,1));
    title("State q(1) : rotary arm angle [rad]");
    legend("Nonlinear", "Bilinear", "EDMD", "Linear");
    savefig("figures\q1.fig");
    saveas(gcf, 'figures\q1.png');
    hold off
    
    figure;
    hold on
    plot(t, comparison.X_nonlinear(2,:,1));
    plot(t, comparison.X_BILIN(2,:,1));
    plot(t, comparison.X_EDMD(2,:,1));
    plot(t, comparison.X_LIN(2,:,1));
    title("State q(2) : pendulum angle [rad]");
    legend("Nonlinear", "Bilinear", "EDMD", "Linear");
    savefig("figures\q2.fig");
    saveas(gcf, 'figures\q2.png');
    hold off
    
    figure;
    hold on
    plot(t, comparison.X_nonlinear(3,:,1));
    plot(t, comparison.X_BILIN(3,:,1));
    plot(t, comparison.X_EDMD(3,:,1));
    plot(t, comparison.X_LIN(3,:,1));
    title("State q(3) : rotary arm speed [rad/s]");
    legend("Nonlinear", "Bilinear", "EDMD", "Linear");
    savefig("figures\q3.fig");
    saveas(gcf, 'figures\q3.png');
    hold off
    
    figure;
    hold on
    plot(t, comparison.X_nonlinear(4,:,1));
    plot(t, comparison.X_BILIN(4,:,1));
    plot(t, comparison.X_EDMD(4,:,1));
    plot(t, comparison.X_LIN(4,:,1));
    title("State q(4) : pendulum speed [rad/s]");
    legend("Nonlinear", "Bilinear", "EDMD", "Linear");
    savefig("figures\q4.fig");
    saveas(gcf, 'figures\q4.png');
    hold off
    
    % 2. Histograms of initial q2 state distribution
    figure;
    histogram(comparison.X_training(2,1,:),BinWidth=2*pi/16,BinLimits=[-pi,pi]);
    title("Training trajectories : distribution of q(2) at t=0 ");
    savefig("figures\q2_distrib_training_trajs.fig");
    saveas(gcf, 'figures\q2_distrib_training_trajs.png');
    
    
    figure;
    histogram(comparison.X_nonlinear(2,1,:),BinWidth=2*pi/16,BinLimits=[-pi,pi]);
    title("Testing trajectories : distribution of q(2) at t=0 ");
    savefig("figures\q2_distrib_testing_trajs.fig");
    saveas(gcf, 'figures\q2_distrib_testing_trajs.png');
    
    
    % All q2 values across all time steps and trajectories
    q2_all_training = reshape(comparison.X_training(2,:,:), 1, []);
    q2_all_testing  = reshape(comparison.X_nonlinear(2,:,:),  1, []);
    
    figure;
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
    y = [comparison.errors.LIN_fro_training comparison.errors.EDMD_fro_training comparison.errors.BILIN_fro_training];
    bar(x,y);
    title("Training errors for original states (Frobenius norm)");
    savefig("figures\training_error_fro.fig");
    saveas(gcf, 'figures\training_error_fro.png');
    
    figure;
    x = ["EDMD" "BILINEAR"];
    y = [comparison.errors.EDMD_fro_training_lifted comparison.errors.BILIN_fro_training_lifted];
    bar(x,y);
    title("Training errors for lifted states (Frobenius norm)");
    savefig("figures\training_error_fro_lifted.fig");
    saveas(gcf, 'figures\training_error_fro_lifted.png');
    
    
    figure;
    x = ["LINEAR" "EDMD" "BILINEAR"];
    y = [comparison.errors.LIN_fro comparison.errors.EDMD_fro comparison.errors.BILIN_fro];
    bar(x,y);
    title("Testing errors for original states (Frobenius norm)");
    savefig("figures\testing_error_original_states_fro.fig");
    saveas(gcf, 'figures\testing_error_original_states_fro.png');
    
    
    figure;
    x = ["EDMD" "BILINEAR"];
    y = [comparison.errors.EDMD_fro_lifted comparison.errors.BILIN_fro_lifted];
    bar(x,y);
    title("Testing errors for lifted states (Frobenius norm)");
    savefig("figures\testing_error_lifted_states_fro.fig");
    saveas(gcf, 'figures\testing_error_lifted_states_fro.png');

    
    figure;
    x = ["LINEAR" "EDMD" "BILINEAR"];
    y = [comparison.errors.LIN_fro_onestep comparison.errors.EDMD_fro_onestep comparison.errors.BILIN_fro_onestep];
    bar(x,y);
    title("Testing errors for original states, one-step prediction (Frobenius norm)");
    savefig("figures\testing_error_original_states_fro_onestep.fig");
    saveas(gcf, 'figures\testing_error_original_states_fro_onestep.png');
    
    
    figure;
    x = ["EDMD" "BILINEAR"];
    y = [comparison.errors.EDMD_fro_lifted_onestep comparison.errors.BILIN_fro_lifted_onestep];
    bar(x,y);
    title("Testing errors for lifted states, one-step prediction (Frobenius norm)");
    savefig("figures\testing_error_lifted_states_fro_onestep.fig");
    saveas(gcf, 'figures\testing_error_lifted_states_fro_onestep.png');
    
    
    
    % 4. Phase portrait of q2 vs q4 (pendulum angle vs pendulum speed)
    q2_reshaped = reshape(comparison.X_training(2,:,:), 1, []);
    q4_reshaped = reshape(comparison.X_training(4,:,:), 1, []);
    figure;
    plot(q2_reshaped, q4_reshaped, '--')
    xlabel('q2 [rad]'); ylabel('q4 [rad/s]');
    title('Phase portrait of q2 vs q4');
    savefig("figures\phase_portrait_q2_vs_q4.fig");
    saveas(gcf, 'figures\phase_portrait_q2_vs_q4.png');
    
    
    % 5. Phase portrait of q1 vs q3 (rotary angle vs rotary speed)
    q1_reshaped = reshape(comparison.X_training(1,:,:), 1, []);
    q3_reshaped = reshape(comparison.X_training(3,:,:), 1, []);
    figure;
    plot(q1_reshaped, q3_reshaped, '--')
    xlabel('q1 [rad]'); ylabel('q3 [rad/s]');
    title('Phase portrait of q1 vs q3');
    savefig("figures\phase_portrait_q1_vs_q3.fig");
    saveas(gcf, 'figures\phase_portrait_q1_vs_q3.png');
    
    
    
    % 6. Input trajectories (to see if suficiently excited)
    figure;
    hold on;
    for i = 1:ntrajs_training
        plot(comparison.U_training(:,:,i));
    end
    title('Training trajectories : input u');
    xlabel('Steps');
    ylabel('u [V]');
    ylim([-16 16]);
    hold off;
    savefig("figures\u_training.fig");
    saveas(gcf, 'figures\u_training.png');
    
    % 7. Singular Values of Phi for the models (to see if some lifting functions are poorly chosen)
    figure;
    semilogy(comparison.svds.LIN);
    title('Singular values of lifted states matrix for linear model');
    savefig("figures\svd_linear.fig");
    saveas(gcf, 'figures\svd_linear.png');
    
    
    figure;
    semilogy(comparison.svds.EDMD);
    title('Singular values of lifted states matrix for EDMD linear model');
    savefig("figures\svd_edmd.fig");
    saveas(gcf, 'figures\svd_edmd.png');
    
    
    figure;
    semilogy(comparison.svds.BILIN);
    title('Singular values of lifted states matrix for EDMD bilinear model');
    savefig("figures\svd_bilinear.fig");
    saveas(gcf, 'figures\svd_bilinear.png');
end