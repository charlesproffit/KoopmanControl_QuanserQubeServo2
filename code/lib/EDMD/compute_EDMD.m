function [EDMD_CT, EDMD_DT] = compute_EDMD(data, f, n, Ts, mode, method)
    if strcmp(mode,'EDMD') || strcmp(mode,'BILINEAR')
        n_used = n.lifted_states;
    else
        n_used = n.states;
    end

    % Step 1 : Compute lifted states
    u = data.training.u;
    q = data.training.q;
    lifted_q_all = zeros(n_used, n.steps, n.trajs_training);
    for i = 1:n.trajs_training
        for j = 1:n.steps
            lifted_q_all(:,j,i) = f(q(:,j,i));
        end
    end


    % Step 2 : Stack data
    U             = reshape(u, n.inputs, []);
    Q             = reshape(q,           n.states, []);
    lifted_Q_all  = reshape(lifted_q_all,    n_used,   []);
    C = Q(:,1:n.steps-1) * pinv(lifted_Q_all(:,1:n.steps-1)); %%WHAT
    normal_idx = [];
    plus_idx  = [];
    for i = 1:n.trajs_training
        base = (i-1)*n.steps;
        normal_idx = [normal_idx, base + (1:n.steps-1)];
        plus_idx  = [plus_idx,  base + (2:n.steps)];
    end

    % Step 3 : Identify models
    if strcmp(mode,'BILINEAR')
        PhiU_all   = lifted_Q_all .* repmat(U, n.lifted_states, 1);
        stacked = [lifted_Q_all(:,normal_idx); PhiU_all(:,normal_idx);]; % removed [lifted_Q; U; PhiU]; % should be equiv to Z' ??
        stacked_phi =  lifted_Q_all(:, normal_idx);
        stacked_plus = [lifted_Q_all(:,plus_idx); PhiU_all(:,plus_idx);]; % should be equiv to Z_Plus ??
        stacked_plus_phi =  lifted_Q_all(:, plus_idx);
        stacked_dot = (stacked_plus - stacked) / Ts;
        stacked_dot_phi  = (stacked_plus_phi - stacked_phi) / Ts;


        % DT identification
        if strcmp(method,'LS')
            ABN_DT = stacked_plus * pinv(stacked);
            stacked_used = stacked_plus;
        elseif strcmp(method,'RIDGE')
            lambda_DT = select_lambda(stacked, stacked_plus_phi, n, [mode '_DT']);
            ABN_DT = stacked_plus_phi * stacked' * (stacked*stacked' + lambda_DT*eye(size(stacked,1)))^(-1);
            EDMD_DT.lambda = lambda_DT;
            stacked_used = stacked_plus_phi;
        end
        A_DT = ABN_DT(1:n_used, 1:n_used);
        N_DT = ABN_DT(1:n_used, n_used+1:end);
        diff_DT       = stacked_used - ABN_DT*stacked;
        diff_DT_phi = diff_DT(1:n_used, :);
        diff_DT       = reshape(diff_DT_phi,       n_used,   n.steps-1, n.trajs_training);
        state_diff_DT = C * diff_DT_phi;
        state_diff_DT = reshape(state_diff_DT, n.states, n.steps-1, n.trajs_training);

        % CT identification
        if strcmp(method,'LS')
            ABN_CT = stacked_dot * pinv(stacked);
            stacked_used = stacked_dot;
        elseif strcmp(method,'RIDGE')
            lambda_CT = select_lambda(stacked, stacked_dot_phi, n, [mode '_CT']);
            ABN_CT = stacked_dot_phi * stacked' * (stacked*stacked' + lambda_CT*eye(size(stacked,1)))^(-1);
            EDMD_CT.lambda = lambda_CT;
            stacked_used = stacked_dot_phi;
        end

        A_CT = ABN_CT(1:n_used, 1:n_used);
        N_CT = ABN_CT(1:n_used, n_used+1:end);
        diff_CT       = stacked_used - ABN_CT*stacked;
        diff_CT_phi = diff_CT(1:n_used, :);
        diff_CT       = reshape(diff_CT_phi,       n_used,   n.steps-1, n.trajs_training);
        state_diff_CT = C * diff_CT_phi;
        state_diff_CT = reshape(state_diff_CT, n.states, n.steps-1, n.trajs_training);

        traj_fro_DT = zeros(n.trajs_training,1);
        traj_fro_CT = zeros(n.trajs_training,1);
        traj_fro_DT_lifted = zeros(n.trajs_training,1);
        traj_fro_CT_lifted = zeros(n.trajs_training,1);
        for i = 1:n.trajs_training
            traj_fro_DT(i)        = norm(state_diff_DT(:,:,i), 'fro');
            traj_fro_CT(i)        = norm(state_diff_CT(:,:,i), 'fro');
            traj_fro_DT_lifted(i) = norm(diff_DT(:,:,i), 'fro');
            traj_fro_CT_lifted(i) = norm(diff_CT(:,:,i), 'fro');
        end

        % Shared fields
        EDMD_DT.condition_number = cond(stacked_phi);
        EDMD_DT.svds = svd(stacked_phi);
        EDMD_CT.condition_number = cond(stacked_phi);
        EDMD_CT.svds = svd(stacked_phi);

        EDMD_DT.A = A_DT;
        EDMD_DT.N = N_DT;
        EDMD_DT.C = C;
        EDMD_CT.A = A_CT;
        EDMD_CT.N = N_CT;
        EDMD_CT.C = C;

        EDMD_DT.training_error_fro        = mean(traj_fro_DT);
        EDMD_DT.training_error_fro_lifted = mean(traj_fro_DT_lifted);
        EDMD_CT.training_error_fro        = mean(traj_fro_CT);
        EDMD_CT.training_error_fro_lifted = mean(traj_fro_CT_lifted);

    else
        stacked = [lifted_Q_all(:,normal_idx); U(:,normal_idx);]; % removed [lifted_Q; U; PhiU]; % should be equiv to Z' ??
        stacked_plus = [lifted_Q_all(:,plus_idx); U(:,plus_idx);]; % should be equiv to Z_Plus ??
        stacked_dot = (stacked_plus - stacked) / Ts;

        % DT identification
        if strcmp(method,'LS')
            AB_DT = stacked_plus * pinv(stacked);
        elseif strcmp(method,'RIDGE')
            lambda_DT = select_lambda(stacked, stacked_plus, n, [mode '_DT']);
            AB_DT = stacked_plus * stacked' * (stacked*stacked' + lambda_DT*eye(size(stacked,1)))^(-1);
            EDMD_DT.lambda = lambda_DT;
        end

        A_DT = AB_DT(1:n_used, 1:n_used);
        B_DT = AB_DT(1:n_used, n_used+1:end);
        EDMD_DT.condition_number = cond(stacked);
        EDMD_DT.svds = svd(stacked);
        diff_DT       = stacked_plus - AB_DT*stacked;
        diff_DT_phi = diff_DT(1:n_used, :);
        diff_DT       = reshape(diff_DT_phi,       n_used,   n.steps-1, n.trajs_training);
        state_diff_DT = C * diff_DT_phi;
        state_diff_DT = reshape(state_diff_DT, n.states, n.steps-1, n.trajs_training);

        % CT identification
        if strcmp(method,'LS')
            AB_CT = stacked_dot * pinv(stacked);
        elseif strcmp(method,'RIDGE')
            lambda_CT = select_lambda(stacked, stacked_dot, n, [mode '_CT']);
            AB_CT = stacked_dot * stacked' * (stacked*stacked' + lambda_CT*eye(size(stacked,1)))^(-1);
            EDMD_CT.lambda = lambda_CT;
        end

        A_CT = AB_CT(1:n_used, 1:n_used);
        B_CT = AB_CT(1:n_used, n_used+1:end);
        EDMD_CT.condition_number = cond(stacked);
        EDMD_CT.svds = svd(stacked);
        diff_CT       = stacked_dot - AB_CT*stacked;
        diff_CT_phi = diff_CT(1:n_used, :);
        diff_CT       = reshape(diff_CT_phi,       n_used,   n.steps-1, n.trajs_training);
        state_diff_CT = C * diff_CT_phi;
        state_diff_CT = reshape(state_diff_CT, n.states, n.steps-1, n.trajs_training);

        traj_fro_DT = zeros(n.trajs_training,1);  traj_fro_CT = zeros(n.trajs_training,1);
        traj_fro_DT_lifted = zeros(n.trajs_training,1);  traj_fro_CT_lifted = zeros(n.trajs_training,1);
        for i = 1:n.trajs_training
            traj_fro_DT(i)        = norm(state_diff_DT(:,:,i), 'fro');
            traj_fro_CT(i)        = norm(state_diff_CT(:,:,i), 'fro');
            traj_fro_DT_lifted(i) = norm(diff_DT(:,:,i), 'fro');
            traj_fro_CT_lifted(i) = norm(diff_CT(:,:,i), 'fro');
        end

        EDMD_DT.A = A_DT;
        EDMD_DT.B = B_DT;
        EDMD_DT.C = C;
        EDMD_CT.A = A_CT;
        EDMD_CT.B = B_CT;
        EDMD_CT.C = C;

        EDMD_DT.training_error_fro        = mean(traj_fro_DT);
        EDMD_DT.training_error_fro_lifted = mean(traj_fro_DT_lifted);
        EDMD_CT.training_error_fro        = mean(traj_fro_CT);
        EDMD_CT.training_error_fro_lifted = mean(traj_fro_CT_lifted);
    end
end

function lambda = select_lambda(stacked, target, n, mode)
    lambdas   = [0, logspace(-4, 4, 50)];
    n_lambdas = length(lambdas);

    traj_fold_assignments = mod(0:n.trajs_training-1, n.cross_val_groups) + 1;
    fold_assignments      = repelem(traj_fold_assignments, n.steps-1);

    val_errors = zeros(n.cross_val_groups, n_lambdas);

    for fold = 1:n.cross_val_groups
        val_mask   = (fold_assignments == fold);
        train_mask = ~val_mask;

        stacked_train = stacked(:, train_mask);
        stacked_val   = stacked(:, val_mask);
        target_train  = target(:, train_mask);
        target_val    = target(:, val_mask);

        for k = 1:n_lambdas
            if lambdas(k) == 0
                AB = target_train * pinv(stacked_train);
            else
                AB = target_train * stacked_train' * (stacked_train*stacked_train' + lambdas(k)*eye(size(stacked_train,1)))^(-1);
            end
            val_errors(fold, k) = norm(target_val - AB*stacked_val, 'fro');
        end
    end

    mean_val_errors = mean(val_errors, 1);
    [~, best_k]     = min(mean_val_errors);
    lambda          = lambdas(best_k);
    fprintf('[%s] Best lambda = %.6f\n', mode, lambda);

    %% Plot
    figure;
    lambda_plot = lambdas;
    lambda_plot(1) = 1e-5;  % replace 0 for log scale display

    semilogx(lambda_plot, mean_val_errors, 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
    xline(max(lambda, 1e-5), 'r--', 'LineWidth', 2, ...
          'DisplayName', sprintf('\\lambda^* = %.2e', lambda));
    xlabel('\lambda');
    ylabel('Mean validation error (Frobenius norm)');
    title(sprintf('Cross-validation curve — %s', mode));
    grid on;
end