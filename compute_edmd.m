function [EDMD, data_lifted] = compute_edmd(data, f, n, mode, method)
    if strcmp(mode,'EDMD')||strcmp(mode,'BILINEAR')
        n_used = n.lifted_states;
    else
        n_used = n.states;
    end
    % Step 1 : Compute lifted states
    lifted_q = zeros(n_used, n.steps+1, n.trajs_training);
    for i = 1:n.trajs_training
        for j = 1:n.steps+1
            lifted_q(:,j,i) = f(data.training.q(:,j,i));
        end
    end
    data_lifted.q = lifted_q;
    data_lifted.u = data.training.u;
    
    % Step 2 : Stack data in the correct format
    q = data.training.q(:,1:n.steps,:);
    lifted_q = data_lifted.q(:,1:n.steps,:);
    lifted_q_plus = data_lifted.q(:,2:n.steps+1,:);
    Q = reshape(q,n.states,[]);
    lifted_Q = reshape(lifted_q,n_used,[]);
    lifted_Q_plus = reshape(lifted_q_plus,n_used,[]);
    U = reshape(data_lifted.u,n.inputs,[]);

    % Step 3 : Compute A, B, C
    if strcmp(mode,'BILINEAR')
        PhiU = lifted_Q.*repmat(U,n.lifted_states,1);
        stacked = [lifted_Q;U;PhiU];
        if strcmp(method,'LS')
            ABN_combined = lifted_Q_plus * pinv(stacked);
        else 
            if strcmp(method,'RIDGE')
                lambda = select_lambda(stacked, lifted_Q_plus, n, mode);
                EDMD.lambda = lambda;
                ABN_combined = lifted_Q_plus * stacked' * (stacked*stacked' + lambda*eye(size(stacked,1)))^(-1);
            else
                disp('Method not implemented');
            end
        end
        A = ABN_combined(:,1:n_used);
        B = ABN_combined(:,n_used+1:n_used+n.inputs);
        N = ABN_combined(:,n_used+n.inputs+1:end);
        C = Q * pinv(lifted_Q);

        % Metrics to identify if the observables are well chosen
        EDMD.condition_number = cond(stacked);
        EDMD.svds = svd(stacked);

        % Error computation
        difference = lifted_Q_plus - A*lifted_Q - B*U - N*PhiU;
        state_diff = C * difference;  % Project residuals to original state space [n.states*(n.steps*n.trajs_training)]
        difference = reshape(difference, n_used, n.steps, n.trajs_training); % [n.states*n.steps*n.trajs_training]
        state_diff = reshape(state_diff, n.states, n.steps, n.trajs_training); % [n.states*n.steps*n.trajs_training]
        % For each trajectory, we compute the Frobenius norm of its error matrix
        traj_fro = zeros(n.trajs_training, 1);
        traj_fro_lifted = zeros(n.trajs_training, 1);
        for i = 1:n.trajs_training
            traj_fro(i) = norm(state_diff(:,:,i), 'fro');
            traj_fro_lifted(i) = norm(difference(:,:,i), 'fro');
        end
        training_error_fro = mean(traj_fro);
        training_error_fro_lifted = mean(traj_fro_lifted);

        EDMD.A = A;
        EDMD.B = B;
        EDMD.C = C;
        EDMD.N = N;
        EDMD.training_error_fro = training_error_fro;
        EDMD.training_error_fro_lifted = training_error_fro_lifted;
    else
        stacked = [lifted_Q;U];
        if strcmp(method,'LS')
            AB_combined = lifted_Q_plus * pinv(stacked);
        else 
            if strcmp(method,'RIDGE')
                lambda = select_lambda(stacked, lifted_Q_plus, n, mode);
                EDMD.lambda = lambda;
                AB_combined = lifted_Q_plus *  stacked' * (stacked*stacked' + lambda*eye(size(stacked,1)))^(-1);
            else
                disp('Method not implemented');
            end
        end
        A = AB_combined(:,1:n_used);
        B = AB_combined(:,n_used+1:end);
        C = Q * pinv(lifted_Q);

        % Metrics to identify if the observables are well chosen
        EDMD.condition_number = cond(stacked);
        EDMD.svds = svd(stacked);        

        % Error computation
        difference = lifted_Q_plus - A*lifted_Q - B*U;
        state_diff = C * difference;  % Project residuals to original state space [n.states*(n.steps*n.trajs_training)]
        difference = reshape(difference, n_used, n.steps, n.trajs_training); % [n.states*n.steps*n.trajs_training]
        state_diff = reshape(state_diff, n.states, n.steps, n.trajs_training); % [n.states*n.steps*n.trajs_training]
        % For each trajectory, we compute the Frobenius norm of its error matrix
        traj_fro = zeros(n.trajs_training, 1);
        traj_fro_lifted = zeros(n.trajs_training, 1);
        for i = 1:n.trajs_training
            traj_fro(i) = norm(state_diff(:,:,i), 'fro');
            traj_fro_lifted(i) = norm(difference(:,:,i), 'fro');
        end
        training_error_fro = mean(traj_fro);
        training_error_fro_lifted = mean(traj_fro_lifted);
    
        EDMD.A = A;
        EDMD.B = B;
        EDMD.C = C;
        EDMD.training_error_fro = training_error_fro;
        EDMD.training_error_fro_lifted = training_error_fro_lifted;
    end

end

function lambda = select_lambda(stacked, lifted_Q_plus, n, mode)

    lambdas = [0, logspace(-4, 4, 50)];
    n_lambdas = length(lambdas);

    traj_fold_assignments = mod(0:n.trajs_training-1, n.cross_val_groups) + 1;
    fold_assignments = repelem(traj_fold_assignments, n.steps);

    val_errors = zeros(n.cross_val_groups, n_lambdas);

    for fold = 1:n.cross_val_groups
        % Split into train / validation
        val_mask   = (fold_assignments == fold);
        train_mask = ~val_mask;

        stacked_train        = stacked(:, train_mask);
        stacked_val          = stacked(:, val_mask);
        lifted_Q_plus_train  = lifted_Q_plus(:, train_mask);
        lifted_Q_plus_val    = lifted_Q_plus(:, val_mask);

        for k = 1:n_lambdas
            if lambdas(k) == 0
                AB = lifted_Q_plus_train * pinv(stacked_train);
            else
                AB = lifted_Q_plus_train * stacked_train' * (stacked_train * stacked_train' + lambdas(k) * eye(size(stacked_train, 1)))^(-1);
            end
            val_errors(fold, k) = norm(lifted_Q_plus_val - AB * stacked_val, 'fro');
        end
    end

    % Average validation error across folds, pick best lambda
    mean_val_errors = mean(val_errors, 1);
    [~, best_k] = min(mean_val_errors);
    lambda = lambdas(best_k);
    fprintf('[%s] Best lambda = %.6f \n', mode, lambda);

    % Plot
    % figure;
    % semilogx(lambdas(2:end), mean_val_errors(2:end), 'o-', 'DisplayName', 'Mean CV error'); hold on;
    % if best_k > 1  % lambda=0 doesn't show on log axis
    %     xline(lambda, 'r--', sprintf('\\lambda = %.4f', lambda));
    % end
    % xlabel('\lambda'); ylabel('Validation error (fro)');
    % title(sprintf('%s: %d-fold crossvalidation for \\lambda selection', mode, n.cross_val_groups));

end