function plots_FL_OL(trajs, Ts, threshold)
    state_names = ["Rotary arm angle [rad]", "Pendulum angle [rad]", ...
                   "Rotary arm speed [rad/s]", "Pendulum speed [rad/s]"];

    n_states = size(trajs.z_nl, 1);
    n_steps  = size(trajs.z_nl, 2);
    n_trajs  = size(trajs.z_nl, 3);
    t        = (0:n_steps-1) * Ts;

    % Identify valid trajectories
    valid = false(1, n_trajs);
    for i = 1:n_trajs
        traj = trajs.z_nl(:,:,i);
        if ~any(isnan(traj(:))) && ~any(isinf(traj(:))) && all(abs(traj(:)) < threshold)
            valid(i) = true;
        end
    end
    fprintf('%d / %d trajectories are valid\n', sum(valid), n_trajs);

    % Input plot
    figure; hold on;
    for i = 1:n_trajs
        if valid(i)
            plot(t, squeeze(trajs.u_nl(1,:,i)));
        end
    end
    xlabel('Time [s]'); ylabel('Input [V]');
    title(sprintf('Inputs — valid trajectories (%d/%d)', sum(valid), n_trajs));
    hold off;

    % State plots: NL vs LTI
    for s = 1:n_states
        figure; hold on;
        h_nl  = NaN; h_lti = NaN;
        for i = 1:n_trajs
            if valid(i)
                h_nl  = plot(t, trajs.z_nl(s,:,i),  'b');
                h_lti = plot(t, trajs.z_lti(s,:,i), 'r--');
            end
        end
        xlabel('Time [s]'); ylabel(state_names(s));
        title(sprintf('State q(%d) — valid trajectories (%d/%d)', s, sum(valid), n_trajs));
        legend([h_nl, h_lti], {'NL', 'LTI'});
        hold off;
    end
end