function plots_FL_CL(trajs, Ts, threshold, gen_title)
    state_names = ["Rotary arm angle [rad]", "Pendulum angle [rad]", ...
                   "Rotary arm speed [rad/s]", "Pendulum speed [rad/s]"];

    n_states = size(trajs.q, 1);
    n_steps  = size(trajs.q, 2);
    n_trajs  = size(trajs.q, 3);
    t        = (0:n_steps-1) * Ts;

    % Identify valid trajectories
    valid = false(1, n_trajs);
    for i = 1:n_trajs
        traj = trajs.q(:,:,i);
        if ~any(isnan(traj(:))) && ~any(isinf(traj(:))) && all(abs(traj(:)) < threshold)
            valid(i) = true;
        end
    end
    fprintf('%d / %d trajectories are valid\n', sum(valid), n_trajs);

    % Input plot
    figure;
    subplot(3,2,1); hold on;
    for i = 1:n_trajs
        if valid(i)
            plot(t, squeeze(trajs.u(1,:,i)));
        end
    end
    title("Input [V]");
    hold off;

    % State plots: NL vs LTI
    for s = 1:n_states
        subplot(3,2,s+1); hold on;
        for i = 1:n_trajs
            if valid(i)
                plot(t, trajs.q(s,:,i),  'b');
            end
        end
        ylabel(state_names(s));
        title(sprintf('q(%d)', s));
        grid on;
    
        if s == n_states
            xlabel('Time [s]');
        end
    
        hold off;
    end
    sgtitle(sprintf('States - %s - valid trajectories (%d/%d)', gen_title, sum(valid), n_trajs));
end