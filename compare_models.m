function comparison = compare_models(data_EDMD, f_lifting, M_BILINEAR_DT, M_BILINEAR_CT, M_EDMD, M_LINEAR, n, Ts)
    q_nonlin = data_EDMD.testing.q;
    q_edmdbilin_DT = zeros(n.lifted_states, n.steps+1, n.trajs_testing);
    q_edmdbilin_CT = zeros(n.lifted_states, n.steps+1, n.trajs_testing);
    q_edmdlin      = zeros(n.lifted_states, n.steps+1, n.trajs_testing);
    q_lin          = zeros(n.states,        n.steps+1, n.trajs_testing);
    u = data_EDMD.testing.u;

    for i = 1:n.trajs_testing
        q_edmdbilin_DT(:,1,i) = f_lifting(q_nonlin(:,1,i));
        q_edmdbilin_CT(:,1,i) = f_lifting(q_nonlin(:,1,i));
        q_edmdlin(:,1,i)      = f_lifting(q_nonlin(:,1,i));
        q_lin(:,1,i)          = q_nonlin(:,1,i);

        for j = 1:n.steps
            ui = u(:,j,i);

            % Bilinear DT
            z = q_edmdbilin_DT(:,j,i);
            q_edmdbilin_DT(:,j+1,i) = M_BILINEAR_DT.A*z + M_BILINEAR_DT.B*ui + M_BILINEAR_DT.N*(z.*ui);

            % Bilinear CT (Euler integration)
            z = q_edmdbilin_CT(:,j,i);
            zdot = M_BILINEAR_CT.A*z + M_BILINEAR_CT.B*ui + M_BILINEAR_CT.N*(z.*ui);
            q_edmdbilin_CT(:,j+1,i) = z + Ts*zdot;

            % Linear EDMD
            q_edmdlin(:,j+1,i) = M_EDMD.A*q_edmdlin(:,j,i) + M_EDMD.B*ui;

            % Linear
            q_lin(:,j+1,i) = M_LINEAR.A*q_lin(:,j,i) + M_LINEAR.B*ui;
        end
    end

    comparison.q_nonlinear   = q_nonlin;
    comparison.q_BILINEAR_DT = q_edmdbilin_DT;
    comparison.q_BILINEAR_CT = q_edmdbilin_CT;
    comparison.q_EDMD        = q_edmdlin;
    comparison.q_LINEAR      = q_lin;
    comparison.u             = u;

    % 1. Testing errors on original states
    err_bilin_DT = zeros(1, n.trajs_testing);
    err_bilin_CT = zeros(1, n.trajs_testing);
    err_edmd     = zeros(1, n.trajs_testing);
    err_lin      = zeros(1, n.trajs_testing);

    for i = 1:n.trajs_testing
        E_bilin_DT = comparison.q_BILINEAR_DT(1:n.states,:,i) - comparison.q_nonlinear(:,:,i);
        E_bilin_CT = comparison.q_BILINEAR_CT(1:n.states,:,i) - comparison.q_nonlinear(:,:,i);
        E_edmd     = comparison.q_EDMD(1:n.states,:,i)        - comparison.q_nonlinear(:,:,i);
        E_lin      = comparison.q_LINEAR(:,:,i)                - comparison.q_nonlinear(:,:,i);
        err_bilin_DT(i) = norm(E_bilin_DT, 'fro');
        err_bilin_CT(i) = norm(E_bilin_CT, 'fro');
        err_edmd(i)     = norm(E_edmd,     'fro');
        err_lin(i)      = norm(E_lin,      'fro');
    end

    comparison.errors.BILINEAR_DT_fro_testing = mean(err_bilin_DT);
    comparison.errors.BILINEAR_CT_fro_testing = mean(err_bilin_CT);
    comparison.errors.EDMD_fro_testing        = mean(err_edmd);
    comparison.errors.LINEAR_fro_testing      = mean(err_lin);

    % 2. Testing errors on lifted states
    err_bilin_DT_lifted = zeros(1, n.trajs_testing);
    err_bilin_CT_lifted = zeros(1, n.trajs_testing);
    err_edmd_lifted     = zeros(1, n.trajs_testing);

    for i = 1:n.trajs_testing
        lifted_true = zeros(n.lifted_states, n.steps+1);
        for t = 1:n.steps+1
            lifted_true(:,t) = f_lifting(comparison.q_nonlinear(:,t,i));
        end
        E_bilin_DT_lifted = comparison.q_BILINEAR_DT(:,:,i) - lifted_true;
        E_bilin_CT_lifted = comparison.q_BILINEAR_CT(:,:,i) - lifted_true;
        E_edmd_lifted     = comparison.q_EDMD(:,:,i)        - lifted_true;
        err_bilin_DT_lifted(i) = norm(E_bilin_DT_lifted, 'fro');
        err_bilin_CT_lifted(i) = norm(E_bilin_CT_lifted, 'fro');
        err_edmd_lifted(i)     = norm(E_edmd_lifted,     'fro');
    end

    comparison.errors.BILINEAR_DT_fro_lifted_testing = mean(err_bilin_DT_lifted);
    comparison.errors.BILINEAR_CT_fro_lifted_testing = mean(err_bilin_CT_lifted);
    comparison.errors.EDMD_fro_lifted_testing        = mean(err_edmd_lifted);

    % 3. Training errors
    comparison.errors.BILINEAR_DT_fro_training        = M_BILINEAR_DT.training_error_fro;
    comparison.errors.BILINEAR_DT_fro_training_lifted = M_BILINEAR_DT.training_error_fro_lifted;
    comparison.errors.BILINEAR_CT_fro_training        = M_BILINEAR_CT.training_error_fro;
    comparison.errors.BILINEAR_CT_fro_training_lifted = M_BILINEAR_CT.training_error_fro_lifted;
    comparison.errors.EDMD_fro_training               = M_EDMD.training_error_fro;
    comparison.errors.EDMD_fro_training_lifted        = M_EDMD.training_error_fro_lifted;
    comparison.errors.LINEAR_fro_training             = M_LINEAR.training_error_fro;
    comparison.errors.LINEAR_fro_training_lifted      = M_LINEAR.training_error_fro_lifted;

    disp('Models compared');
end