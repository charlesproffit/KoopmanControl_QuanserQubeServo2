function comparison = compare_FL(trajs_cell, labels, Ts, threshold, varargin)

% Optional : LQR matrices used for all trajs (the same for all models)
%   Q
%   R

    compute_lqr = false;
    n_metrics = 3;
    if nargin > 4
        n_metrics = 4;
        compute_lqr = true;
        Q = varargin{1};
        R = varargin{2};
    end

    n_models = length(trajs_cell);
    n_steps  = size(trajs_cell{1}.X, 2);
    n_states = size(trajs_cell{1}.X, 1);
    n_trajs  = size(trajs_cell{1}.X, 3);
    t        = (0:n_steps-1) * Ts;

    colors = lines(n_models);

    %% 1. Validity mask per model
    valid = false(n_models, n_trajs);
    for m = 1:n_models
        for i = 1:n_trajs
            traj = trajs_cell{m}.X(:,:,i);
            if ~any(isnan(traj(:))) && ~any(isinf(traj(:))) && all(abs(traj(:)) < threshold)
                valid(m,i) = true;
            end
        end
        fprintf('[%s] %d / %d valid trajectories\n', labels{m}, sum(valid(m,:)), n_trajs);
    end

    %% 2. Integral state error per model
    errors = zeros(n_models, n_trajs);
    for m = 1:n_models
        for i = 1:n_trajs
            if valid(m,i)
                e = trajs_cell{m}.X(:,:,i) - trajs_cell{m}.Ref(:,:,i);
                errors(m,i) = Ts * sum(vecnorm(e, 2, 1));
            end
        end
    end
    comparison.errors = errors;
    comparison.labels = labels;

    mean_errors = zeros(n_models, 1);
    for m = 1:n_models
        valid_errors = errors(m, valid(m,:));
        mean_errors(m) = mean(valid_errors);
        fprintf('[%s] Mean integral state error = %.4f\n', labels{m}, mean_errors(m));
    end
    comparison.mean_errors = mean_errors;

    %% 3. Comparison plot
    ref_traj_idx = find(valid(1,:), 1);
    state_names = ["q_1 [rad]", "q_2 [rad]", "q_3 [rad/s]", "q_4 [rad/s]"];

    figure('Name', 'CL Comparison : single trajectory');
    for s = 1:n_states
        subplot(n_states+1, 1, s); hold on;
        for m = 1:n_models
            if valid(m, ref_traj_idx)
                plot(t, trajs_cell{m}.X(s,:,ref_traj_idx), 'Color', colors(m,:), 'LineWidth', 2);
            end
        end
        plot(t, trajs_cell{1}.Ref(s,:,ref_traj_idx), 'k--', 'LineWidth', 1);
        ylabel(state_names(s));
        grid on; hold off;
        if s == 1
            legend([labels, {'ref'}], 'Location', 'best');
        end
    end
    subplot(n_states+1, 1, n_states+1); hold on;
    for m = 1:n_models
        if valid(m, ref_traj_idx)
            plot(t, trajs_cell{m}.U(1,:,ref_traj_idx), 'Color', colors(m,:), 'LineWidth', 2);
        end
    end
    ylabel('u [V]'); xlabel('Time [s]');
    grid on; hold off;
    sgtitle('Closed-loop comparison');

    %% 4. Individual plots
    for m = 1:n_models
        figure('Name', sprintf('CL : %s', labels{m}));
        for s = 1:n_states
            subplot(n_states+1, 1, s); hold on;
            for i = 1:n_trajs
                if valid(m,i)
                    plot(t, trajs_cell{m}.X(s,:,i), 'LineWidth', 1);
                end
            end
            plot(t, trajs_cell{m}.Ref(s, :, find(valid(m,:),1)), 'k--', 'LineWidth', 1);
            ylabel(state_names(s));
            grid on; hold off;
        end
        subplot(n_states+1, 1, n_states+1); hold on;
        for i = 1:n_trajs
            if valid(m,i)
                plot(t, trajs_cell{m}.U(1,:,i), 'LineWidth', 1);
            end
        end
        ylabel('u [V]'); xlabel('Time [s]');
        grid on; hold off;
        sgtitle(sprintf('%s : %d/%d valid, mean error = %.4f', ...
            labels{m}, sum(valid(m,:)), n_trajs, mean_errors(m)));
    end

    % %% 5a. Error bar chart
    % figure('Name', 'Mean integral state error');
    % bar(mean_errors);
    % set(gca, 'XTickLabel', labels);
    % ylabel('Mean integral state error');
    % title('CL performance comparison');
    % grid on;

    %% 5b. Control effort per model
    effort = zeros(n_models, n_trajs);
    for m = 1:n_models
        for i = 1:n_trajs
            if valid(m,i)
                effort(m,i) = Ts * sum(abs(trajs_cell{m}.U(1,:,i)));
            end
        end
    end
    comparison.effort = effort;
    
    mean_effort = zeros(n_models, 1);
    for m = 1:n_models
        valid_effort = effort(m, valid(m,:));
        mean_effort(m) = mean(valid_effort);
        fprintf('[%s] Mean control effort = %.4f\n', labels{m}, mean_effort(m));
    end
    comparison.mean_effort = mean_effort;
    
    %% 5c. Controlled states error (q1 and q3 only)
    controlled_idx = [1, 3];
    errors_controlled = zeros(n_models, n_trajs);
    for m = 1:n_models
        for i = 1:n_trajs
            if valid(m,i)
                e = trajs_cell{m}.X(controlled_idx,:,i) - trajs_cell{m}.Ref(controlled_idx,:,i);
                errors_controlled(m,i) = Ts * sum(vecnorm(e, 2, 1));
            end
        end
    end
    comparison.errors_controlled = errors_controlled;
    
    mean_errors_controlled = zeros(n_models, 1);
    for m = 1:n_models
        valid_ec = errors_controlled(m, valid(m,:));
        mean_errors_controlled(m) = mean(valid_ec);
        fprintf('[%s] Mean controlled state error (q1,q3) = %.4f\n', labels{m}, mean_errors_controlled(m));
    end
    comparison.mean_errors_controlled = mean_errors_controlled;

    %% 6. LQR objective (optional)
    if compute_lqr
        lqr_costs = zeros(n_models, n_trajs);
        for m = 1:n_models
            for i = 1:n_trajs
                if valid(m,i)
                    J_traj = 0;
                    for k = 1:n_steps
                        x_k = trajs_cell{m}.X(:,k,i);
                        r_k = trajs_cell{m}.Ref(:,k,i);
                        u_k = trajs_cell{m}.U(:,k,i);
                        e_k = x_k - r_k;
                        J_traj = J_traj + e_k'*Q*e_k + u_k'*R*u_k;
                    end
                    lqr_costs(m,i) = J_traj;
                end
            end
        end
        comparison.lqr_costs = lqr_costs;

        mean_lqr = zeros(n_models, 1);
        for m = 1:n_models
            valid_costs = lqr_costs(m, valid(m,:));
            mean_lqr(m) = mean(valid_costs);
            fprintf('[%s] Mean LQR cost = %.4f\n', labels{m}, mean_lqr(m));
        end
        comparison.mean_lqr = mean_lqr;
    end

    %% 6. Bar charts
    figure('Name', 'CL performance summary');
    
    subplot(1,n_metrics,1);
    bar(mean_errors);
    set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 15);
    ylabel('Mean integral state error (all states)');
    title('Full state error');
    grid on;
    
    subplot(1,n_metrics,2);
    bar(mean_errors_controlled);
    set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 15);
    ylabel('Mean integral state error (q1, q3)');
    title('Controlled state error');
    grid on;
    
    subplot(1,n_metrics,3);
    bar(mean_effort);
    set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 15);
    ylabel('Mean integral control effort');
    title('Control effort');
    grid on;

    if compute_lqr
        subplot(1, n_metrics, 4);
        bar(mean_lqr);
        set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 15);
        ylabel('Mean LQR cost J');
        title('LQR objective');
        grid on;
    end
    
    sgtitle('Closed-loop performance comparison');
end