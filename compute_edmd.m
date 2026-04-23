function [EDMD_CT, EDMD_DT, data_lifted] = compute_edmd(data, f, n, Ts, mode, method)
    if strcmp(mode,'EDMD') || strcmp(mode,'BILINEAR')
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

    % Step 2 : Stack data
    q          = data.training.q(:,1:n.steps,:);
    lifted_q   = data_lifted.q(:,1:n.steps,:);
    lifted_q_plus = data_lifted.q(:,2:n.steps+1,:);
    lifted_q_dot  = (lifted_q_plus - lifted_q) / Ts;

    Q             = reshape(q,           n.states, []);
    lifted_Q      = reshape(lifted_q,    n_used,   []);
    lifted_Q_plus = reshape(lifted_q_plus, n_used, []);
    lifted_Q_dot  = reshape(lifted_q_dot,  n_used, []);
    U             = reshape(data_lifted.u, n.inputs, []);
    C = Q * pinv(lifted_Q);

    % Step 3 : Identify models
    if strcmp(mode,'BILINEAR')
        PhiU   = lifted_Q .* repmat(U, n.lifted_states, 1);
        stacked = [lifted_Q; U; PhiU];

        % DT identification
        if strcmp(method,'LS')
            ABN_DT = lifted_Q_plus * pinv(stacked);
        elseif strcmp(method,'RIDGE')
            lambda_DT = select_lambda(stacked, lifted_Q_plus, n, [mode '_DT']);
            ABN_DT = lifted_Q_plus * stacked' * (stacked*stacked' + lambda_DT*eye(size(stacked,1)))^(-1);
            EDMD_DT.lambda = lambda_DT;
        end

        A_DT = ABN_DT(:, 1:n_used);
        B_DT = ABN_DT(:, n_used+1:n_used+n.inputs);
        N_DT = ABN_DT(:, n_used+n.inputs+1:end);
        diff_DT       = lifted_Q_plus - A_DT*lifted_Q - B_DT*U - N_DT*PhiU;
        state_diff_DT = C * diff_DT;
        diff_DT       = reshape(diff_DT,       n_used,   n.steps, n.trajs_training);
        state_diff_DT = reshape(state_diff_DT, n.states, n.steps, n.trajs_training);

        % CT identification
        if strcmp(method,'LS')
            ABN_CT = lifted_Q_dot * pinv(stacked);
        elseif strcmp(method,'RIDGE')
            lambda_CT = select_lambda(stacked, lifted_Q_dot, n, [mode '_CT']);
            ABN_CT = lifted_Q_dot * stacked' * (stacked*stacked' + lambda_CT*eye(size(stacked,1)))^(-1);
            EDMD_CT.lambda = lambda_CT;
        end

        A_CT = ABN_CT(:, 1:n_used);
        B_CT = ABN_CT(:, n_used+1:n_used+n.inputs);
        N_CT = ABN_CT(:, n_used+n.inputs+1:end);
        diff_CT       = lifted_Q_dot - A_CT*lifted_Q - B_CT*U - N_CT*PhiU;
        state_diff_CT = C * diff_CT;
        diff_CT       = reshape(diff_CT,       n_used,   n.steps, n.trajs_training);
        state_diff_CT = reshape(state_diff_CT, n.states, n.steps, n.trajs_training);

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
        EDMD_DT.condition_number = cond(stacked);
        EDMD_DT.svds = svd(stacked);
        EDMD_CT.condition_number = cond(stacked);
        EDMD_CT.svds = svd(stacked);

        EDMD_DT.A = A_DT;
        EDMD_DT.B = B_DT;
        EDMD_DT.N = N_DT;
        EDMD_DT.C = C;
        EDMD_CT.A = A_CT;
        EDMD_CT.B = B_CT;
        EDMD_CT.N = N_CT;
        EDMD_CT.C = C;

        EDMD_DT.training_error_fro        = mean(traj_fro_DT);
        EDMD_DT.training_error_fro_lifted = mean(traj_fro_DT_lifted);
        EDMD_CT.training_error_fro        = mean(traj_fro_CT);
        EDMD_CT.training_error_fro_lifted = mean(traj_fro_CT_lifted);

    else
        stacked = [lifted_Q; U];

        % DT identification
        if strcmp(method,'LS')
            AB_DT = lifted_Q_plus * pinv(stacked);
        elseif strcmp(method,'RIDGE')
            lambda_DT = select_lambda(stacked, lifted_Q_plus, n, [mode '_DT']);
            AB_DT = lifted_Q_plus * stacked' * (stacked*stacked' + lambda_DT*eye(size(stacked,1)))^(-1);
            EDMD_DT.lambda = lambda_DT;
        end

        A_DT = AB_DT(:, 1:n_used);
        B_DT = AB_DT(:, n_used+1:end);
        EDMD_DT.condition_number = cond(stacked);
        EDMD_DT.svds = svd(stacked);
        diff_DT       = lifted_Q_plus - A_DT*lifted_Q - B_DT*U;
        state_diff_DT = C * diff_DT;
        diff_DT       = reshape(diff_DT,       n_used,   n.steps, n.trajs_training);
        state_diff_DT = reshape(state_diff_DT, n.states, n.steps, n.trajs_training);

        % CT identification
        if strcmp(method,'LS')
            AB_CT = lifted_Q_dot * pinv(stacked);
        elseif strcmp(method,'RIDGE')
            lambda_CT = select_lambda(stacked, lifted_Q_dot, n, [mode '_CT']);
            AB_CT = lifted_Q_dot * stacked' * (stacked*stacked' + lambda_CT*eye(size(stacked,1)))^(-1);
            EDMD_CT.lambda = lambda_CT;
        end

        A_CT = AB_CT(:, 1:n_used);
        B_CT = AB_CT(:, n_used+1:end);
        EDMD_CT.condition_number = cond(stacked);
        EDMD_CT.svds = svd(stacked);
        diff_CT       = lifted_Q_dot - A_CT*lifted_Q - B_CT*U;
        state_diff_CT = C * diff_CT;
        diff_CT       = reshape(diff_CT,       n_used,   n.steps, n.trajs_training);
        state_diff_CT = reshape(state_diff_CT, n.states, n.steps, n.trajs_training);

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
    fold_assignments      = repelem(traj_fold_assignments, n.steps);

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
                AB = target_train * stacked_train' * ...
                     (stacked_train*stacked_train' + lambdas(k)*eye(size(stacked_train,1)))^(-1);
            end
            val_errors(fold, k) = norm(target_val - AB*stacked_val, 'fro');
        end
    end

    mean_val_errors = mean(val_errors, 1);
    [~, best_k]     = min(mean_val_errors);
    lambda          = lambdas(best_k);
    fprintf('[%s] Best lambda = %.6f\n', mode, lambda);
end