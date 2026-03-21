function [EDMD, data_lifted] = compute_edmd(data, f, n, mode)
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
    % Metrics to identify if the state space is well covered
    EDMD.condition_number = cond(lifted_Q);
    EDMD.svds = svd(lifted_Q);

    % Step 3 : Compute A, B, C
    if strcmp(mode,'BILINEAR')
        PhiU = lifted_Q.*repmat(U,n.lifted_states,1);
        ABN_combined = lifted_Q_plus * pinv([lifted_Q;U;PhiU]);
        A = ABN_combined(:,1:n_used);
        B = ABN_combined(:,n_used+1:n_used+n.inputs);
        N = ABN_combined(:,n_used+n.inputs+1:end);
        C = Q * pinv(lifted_Q);

        % Error computation
        difference = lifted_Q_plus - A*lifted_Q - B*U - N*PhiU;
        state_diff = C * difference;  % Project residuals to original state space [n.states*(n.steps*n.trajs_training)]
        state_diff = reshape(state_diff, n.states, n.steps, n.trajs_training); % [n.states*n.steps*n.trajs_training]
        % For each trajectory, we compute the Frobenius norm of its error matrix
        traj_fro = zeros(n.trajs_training, 1);
        for i = 1:n.trajs_training
            traj_fro(i) = norm(state_diff(:,:,i), 'fro');
        end
        training_error_fro = mean(traj_fro); % Then we average
        % training_error_2norm = norm(difference, 2);
        % training_error_fro = norm(difference, 'fro');

        EDMD.A = A;
        EDMD.B = B;
        EDMD.C = C;
        EDMD.N = N;
        % EDMD.training_error_2norm = training_error_2norm;
        EDMD.training_error_fro = training_error_fro;
    else
        AB_combined = lifted_Q_plus * pinv([lifted_Q;U]);
        A = AB_combined(:,1:n_used);
        B = AB_combined(:,n_used+1:end);
        C = Q * pinv(lifted_Q);

        % Error computation
        difference = lifted_Q_plus - A*lifted_Q - B*U;
        state_diff = C * difference;  % Project residuals to original state space [n.states*(n.steps*n.trajs_training)]
        state_diff = reshape(state_diff, n.states, n.steps, n.trajs_training); % [n.states*n.steps*n.trajs_training]
        % For each trajectory, we compute the Frobenius norm of its error matrix
        traj_fro = zeros(n.trajs_training, 1);
        for i = 1:n.trajs_training
            traj_fro(i) = norm(state_diff(:,:,i), 'fro');
        end
        training_error_fro = mean(traj_fro); % Then we average
        % training_error_2norm = norm(difference, 2); % vecnorm
        % training_error_fro = norm(difference, 'fro');
    
        EDMD.A = A;
        EDMD.B = B;
        EDMD.C = C;
        % EDMD.training_error_2norm = training_error_2norm;
        EDMD.training_error_fro = training_error_fro;
    end

end