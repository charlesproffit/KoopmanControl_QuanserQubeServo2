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

    % Step 3 : Compute A, B, C
    if strcmp(mode,'BILINEAR')
        ABN_combined = lifted_Q_plus * pinv([lifted_Q;U;lifted_Q.*repmat(U,n.lifted_states,1)]);
        A = ABN_combined(:,1:n_used);
        B = ABN_combined(:,n_used+1:n_used+n.inputs);
        N = ABN_combined(:,n_used+n.inputs+1:end);
        C = Q * pinv(lifted_Q);
        difference = lifted_Q_plus - A*lifted_Q - B*U - N*lifted_Q.*repmat(U,n.lifted_states,1);
        training_error_2norm = norm(difference, 2);
        training_error_fro = norm(difference, 'fro');

        EDMD.A = A;
        EDMD.B = B;
        EDMD.C = C;
        EDMD.N = N;
        EDMD.training_error_2norm = training_error_2norm;
        EDMD.training_error_fro = training_error_fro;
    else
        AB_combined = lifted_Q_plus * pinv([lifted_Q;U]);
        A = AB_combined(:,1:n_used);
        B = AB_combined(:,n_used+1:end);
        C = Q * pinv(lifted_Q);
        difference = lifted_Q_plus - A*lifted_Q - B*U;
        training_error_2norm = norm(difference, 2);
        training_error_fro = norm(difference, 'fro');
    
        EDMD.A = A;
        EDMD.B = B;
        EDMD.C = C;
        EDMD.training_error_2norm = training_error_2norm;
        EDMD.training_error_fro = training_error_fro;
    end

end