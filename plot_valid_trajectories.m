function plot_valid_trajectories(trajs, Ts, threshold)
    n_states = size(trajs.q, 1);
    n_steps  = size(trajs.q, 2);
    n_trajs  = size(trajs.q, 3);
    t = (0:n_steps-1) * Ts;

    state_names = ["Rotary arm angle [rad]","Pendulum angle [rad]", "Rotary arm speed [rad/s]", "Pendulum speed [rad/s]"];

    % Identify valid trajectories
    valid = false(1, n_trajs);
    for i = 1:n_trajs
        traj = trajs.q(:,:,i);
        
        if ~any(isnan(traj(:))) && ~any(isinf(traj(:))) && all(abs(traj(:)) < threshold)
            valid(i) = true;
        end
    end

    fprintf('%d / %d trajectories are valid\n', sum(valid), n_trajs);

    for s = 1:n_states
        figure; hold on;
        for i = 1:n_trajs
            if valid(i)
                plot(t, trajs.q(s,:,i));
            end
        end
        xlabel('Time [s]');
        ylabel(state_names(s));
        title(sprintf('State q(%d) — valid trajectories (%d/%d)', s, sum(valid), n_trajs));
        hold off;
    end
end