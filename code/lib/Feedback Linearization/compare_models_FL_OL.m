function trajs = compare_models_FL_OL(M_MBFL, M_DDFL, f_discrete, n, Ts)
    %% 1. Sweep over each state
    state_names = {'q_1', 'q_2', 'q_3', 'q_4'};
    state_ranges = {linspace(-pi,   pi,   300), linspace(-pi,   pi,   300), linspace(-5*pi, 5*pi, 300), linspace(-5*pi, 5*pi, 300)};

    q_base = [0; 0; 0; 0];
    gamma_ddfl_sweeps = zeros(n.states, 300);
    gamma_mbfl_sweeps = zeros(n.states, 300);
    etta_ddfl_sweeps  = zeros(n.states, 300);
    etta_mbfl_sweeps  = zeros(n.states, 300);

    for s = 1:n.states
        vec = state_ranges{s};
        for k = 1:length(vec)
            q = q_base;
            q(s) = vec(k);
            gamma_ddfl_sweeps(s,k) = M_DDFL.gamma(q);
            gamma_mbfl_sweeps(s,k) = M_MBFL.gamma(q);
            etta_ddfl_sweeps(s,k)  = M_DDFL.etta(q);
            etta_mbfl_sweeps(s,k)  = M_MBFL.etta(q);
        end
    end

    %% 2. Gamma comparison
    figure('Name', 'Gamma comparison: MBFL vs DDFL');
    for s = 1:n.states
        subplot(2, 2, s);
        vec = state_ranges{s};
        plot(vec, gamma_mbfl_sweeps(s,:), 'b-',  'LineWidth', 2); hold on;
        plot(vec, gamma_ddfl_sweeps(s,:), 'r--', 'LineWidth', 2);
        xlabel(state_names{s}); ylabel('\gamma');
        title(sprintf('Sweep over %s', state_names{s}));
        legend('Model-based', 'DDFL');
        grid on; hold off;
    end
    sgtitle('\gamma comparison: MBFL vs DDFL (other states = 0)');

    %% 3. Etta comparison
    figure('Name', 'Etta comparison: MBFL vs DDFL');
    for s = 1:n.states
        subplot(2, 2, s);
        vec = state_ranges{s};
        plot(vec, etta_mbfl_sweeps(s,:), 'b-',  'LineWidth', 2); hold on;
        plot(vec, etta_ddfl_sweeps(s,:), 'r--', 'LineWidth', 2);
        xlabel(state_names{s}); ylabel('\eta');
        title(sprintf('Sweep over %s', state_names{s}));
        legend('Model-based', 'DDFL');
        grid on; hold off;
    end
    sgtitle('\eta comparison: MBFL vs DDFL (other states = 0)');

    %% 4. Comparison of LTI Linearized System with original NL system
    noise = 2*n.umax*rand(n.inputs, n.steps, n.trajs_control) - n.umax;

    q_nl  = zeros(n.states, n.steps, n.trajs_control);
    z_nl_ddfl  = zeros(n.controlled_states, n.steps, n.trajs_control);
    z_nl_mbfl  = zeros(n.controlled_states, n.steps, n.trajs_control);
    z_lti_ddfl = zeros(n.controlled_states, n.steps, n.trajs_control);
    z_lti_mbfl = zeros(n.controlled_states, n.steps, n.trajs_control);
    u_nl  = zeros(n.inputs, n.steps, n.trajs_control);

    for i = 1:n.trajs_control
        % q0
        q_nl(1,1,i) = 0;
        q_nl(2,1,i) = 2*(2*rand() - 1);
        q_nl(3,1,i) = 10*(2*rand() - 1);
        q_nl(4,1,i) = 10*(2*rand() - 1);

        % z0 = T(q0)
        z_nl_ddfl(:,1,i)  = M_DDFL.T(q_nl(:,1,i));
        z_nl_mbfl(:,1,i)  = M_MBFL.T(q_nl(:,1,i));
        z_lti_ddfl(:,1,i) = M_DDFL.T(q_nl(:,1,i));
        z_lti_mbfl(:,1,i) = M_MBFL.T(q_nl(:,1,i));


        for j = 1:n.steps-1
            % random u applied to NL system => simulate z_nl traj
            u = noise(:,j,i);
            u_nl(:,j,i)      = u;
            q_nl(:,j+1,i)    = f_discrete(q_nl(:,j,i), u);

            % lift NL trajectory
            z_nl_ddfl(:,j+1,i) = M_DDFL.T(q_nl(:,j+1,i));
            z_nl_mbfl(:,j+1,i) = M_MBFL.T(q_nl(:,j+1,i));

            % for ddfl
            % v = (u - etta(q)) / gamma(q) and simulate z_lti traj
            ett  = M_DDFL.etta(q_nl(:,j,i));
            gamm = M_DDFL.gamma(q_nl(:,j,i));
            v    = (u - ett) / gamm;
            z_lti_dot = M_DDFL.A * z_lti_ddfl(:,j,i) + M_DDFL.B * v;
            z_lti_ddfl(:,j+1,i) = z_lti_ddfl(:,j,i)+Ts*z_lti_dot;

            % for mbfl 
            % v = (u - etta(q)) / gamma(q) and simulate z_lti traj
            ett  = M_MBFL.etta(q_nl(:,j,i));
            gamm = M_MBFL.gamma(q_nl(:,j,i));
            v    = (u - ett) / gamm;
            z_lti_dot = M_MBFL.A * z_lti_mbfl(:,j,i) + M_MBFL.B * v;
            z_lti_mbfl(:,j+1,i) = z_lti_mbfl(:,j,i)+Ts*z_lti_dot;
        end
    end

    trajs.q_nl  = q_nl;
    trajs.z_nl_mbfl  = z_nl_mbfl;
    trajs.z_nl_ddfl  = z_nl_ddfl;
    trajs.z_lti_mbfl = z_lti_mbfl;
    trajs.z_lti_ddfl = z_lti_ddfl;
    trajs.u_nl  = u_nl;
    disp('Trajectories for comparison collected');


    %% 5. Plots 
    threshold = 1000;

    % z-space state names (T maps to 2D space)
    z_names = ["z_1 (T_1(q))", "z_2 (T_2(q))"];

    t = (0:n.steps-1) * Ts;

    % Identify valid trajectories
    valid = false(1, n.trajs_control);
    for i = 1:n.trajs_control
        traj_nl_ddfl  = trajs.z_nl_ddfl(:,:,i);
        traj_nl_mbfl  = trajs.z_nl_mbfl(:,:,i);
        traj_lti_mbfl = trajs.z_lti_mbfl(:,:,i);
        traj_lti_ddfl = trajs.z_lti_ddfl(:,:,i);

        if ~any(isnan(traj_nl_ddfl(:)))  && ~any(isinf(traj_nl_ddfl(:)))  && all(abs(traj_nl_ddfl(:))  < threshold) && ...
           ~any(isnan(traj_nl_mbfl(:)))  && ~any(isinf(traj_nl_mbfl(:)))  && all(abs(traj_nl_mbfl(:))  < threshold) && ...
           ~any(isnan(traj_lti_mbfl(:))) && ~any(isinf(traj_lti_mbfl(:))) && all(abs(traj_lti_mbfl(:)) < threshold) && ...
           ~any(isnan(traj_lti_ddfl(:))) && ~any(isinf(traj_lti_ddfl(:))) && all(abs(traj_lti_ddfl(:)) < threshold)
            valid(i) = true;
        end
    end
    fprintf('%d / %d trajectories are valid\n', sum(valid), n.trajs_control);

    % z_nl vs z_lti per component
    figure;
    for s = 1:n.controlled_states
        subplot(n.controlled_states,1,s);
        hold on;
        h_nl_mbfl  = NaN;
        h_nl_ddfl  = NaN;
        h_lti_mbfl = NaN;
        h_lti_ddfl = NaN;
        for i = 1:1%n.trajs_control
            if valid(i)
                h_nl_mbfl  = plot(t, trajs.z_nl_mbfl(s,:,i),  'b',   'LineWidth', 0.8);
                h_nl_ddfl  = plot(t, trajs.z_nl_ddfl(s,:,i),  'r',   'LineWidth', 0.8);
                h_lti_mbfl = plot(t, trajs.z_lti_mbfl(s,:,i), 'b--', 'LineWidth', 0.8);
                h_lti_ddfl = plot(t, trajs.z_lti_ddfl(s,:,i), 'r--', 'LineWidth', 0.8);
            end
        end
        xlabel('Time [s]');
        ylabel(z_names(s));
        title(sprintf('%s — NL lifted vs LTI (%d/%d valid)', z_names(s), sum(valid), n.trajs_control));
        legend([h_nl_mbfl, h_nl_ddfl, h_lti_mbfl, h_lti_ddfl], {'z_{NL, MBFL} = T(q)', 'z_{NL, DDFL} = T(q)', 'z_{LTI, MBFL}', 'z_{LTI, DDFL}'});
        grid on; hold off;
    end
end