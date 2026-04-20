function comparison = compare_models(data_EDMD, f_lifting, M_BILINEAR, M_EDMD, M_LINEAR, n)
    q_nonlin = data_EDMD.testing.q;
    q_edmdbilin = zeros(n.lifted_states, n.steps+1, n.trajs_testing);
    q_edmdlin = zeros(n.lifted_states, n.steps+1, n.trajs_testing);
    q_lin = zeros(n.states, n.steps+1, n.trajs_testing);
    u = data_EDMD.testing.u;

    % Reminder : States q :
    %   q(1) => Angle of rotary arm
    %   q(2) => Angle of pendulum arm
    %   q(3) => Angular speed of rotary arm q(1)_dot
    %   q(4) => Angular speed of pendulum arm q(2)_dot

    for i=1:n.trajs_testing
        % Initial state shared between trajs and models
        q_edmdbilin(:,1,i) = f_lifting(q_nonlin(:,1,i)); % Lifted EDMD Bilinear
        q_edmdlin(:,1,i) = f_lifting(q_nonlin(:,1,i)); % Lifted EDMD Linear
        q_lin(:,1,i) = q_nonlin(:,1,i); % Not lifted Linear

        for j = 1:n.steps
            % Bilinear EDMD model (lifting, computing evolution then going back to original states)
            q_edmdbilin(:,j+1,i) = M_BILINEAR.A*q_edmdbilin(:,j,i) + M_BILINEAR.B*u(:,j,i) + M_BILINEAR.N*q_edmdbilin(:,j,i)*u(:,j,i);
            
            % Linear EDMD model (lifting, computing evolution then going back to original states)
            q_edmdlin(:,j+1,i) = M_EDMD.A*q_edmdlin(:,j,i) + M_EDMD.B*u(:,j,i);
            
            % Linear model
            q_lin(:,j+1,i) = M_LINEAR.A*q_lin(:,j,i) + M_LINEAR.B*u(:,j,i);
        end
    end
    comparison.q_nonlinear = q_nonlin(:,:,:);
    comparison.q_BILINEAR = q_edmdbilin(:,:,:);
    comparison.q_EDMD = q_edmdlin(:,:,:);
    comparison.q_LINEAR = q_lin(:,:,:);
    comparison.u = u;
    
    % I compute various errors (don't know which one is better)    
    % 1. Testing error on original states
    err_bilin = zeros(1, n.trajs_testing);
    err_edmd  = zeros(1, n.trajs_testing);
    err_lin   = zeros(1, n.trajs_testing);
    
    for i = 1:n.trajs_testing
        E_bilin = comparison.q_BILINEAR(1:n.states,:,i) - comparison.q_nonlinear(:,:,i);
        E_edmd  = comparison.q_EDMD(1:n.states,:,i)     - comparison.q_nonlinear(:,:,i);
        E_lin   = comparison.q_LINEAR(:,:,i)             - comparison.q_nonlinear(:,:,i);
        err_bilin(i) = norm(E_bilin, 'fro');
        err_edmd(i)  = norm(E_edmd,  'fro');
        err_lin(i)   = norm(E_lin,   'fro');
    end
    
    comparison.errors.BILINEAR_fro_testing = mean(err_bilin);
    comparison.errors.EDMD_fro_testing     = mean(err_edmd);
    comparison.errors.LINEAR_fro_testing   = mean(err_lin);
    
    % 2. Testing errors on lifted states
    err_bilin_lifted = zeros(1, n.trajs_testing);
    err_edmd_lifted  = zeros(1, n.trajs_testing);
    
    for i = 1:n.trajs_testing
        lifted_true = zeros(n.lifted_states, n.steps+1);
        for t = 1:n.steps+1
            lifted_true(:,t) = f_lifting(comparison.q_nonlinear(:,t,i));
        end
        E_bilin_lifted = comparison.q_BILINEAR(:,:,i) - lifted_true;
        E_edmd_lifted  = comparison.q_EDMD(:,:,i)     - lifted_true;
        err_bilin_lifted(i) = norm(E_bilin_lifted, 'fro');
        err_edmd_lifted(i)  = norm(E_edmd_lifted,  'fro');
    end
    
    comparison.errors.BILINEAR_fro_lifted_testing = mean(err_bilin_lifted);
    comparison.errors.EDMD_fro_lifted_testing     = mean(err_edmd_lifted);

    %4. Training erros
    comparison.errors.BILINEAR_fro_training_lifted = M_BILINEAR.training_error_fro_lifted;
    comparison.errors.BILINEAR_fro_training = M_BILINEAR.training_error_fro;
    comparison.errors.EDMD_fro_training_lifted = M_EDMD.training_error_fro_lifted;
    comparison.errors.EDMD_fro_training = M_EDMD.training_error_fro;
    comparison.errors.LINEAR_fro_training_lifted = M_LINEAR.training_error_fro_lifted;
    comparison.errors.LINEAR_fro_training = M_LINEAR.training_error_fro;

    disp('Models compared');
end