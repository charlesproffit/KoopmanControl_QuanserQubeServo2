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
        step_errors_bilin = zeros(n.steps+1, 1);
        step_errors_edmd  = zeros(n.steps+1, 1);
        step_errors_lin   = zeros(n.steps+1, 1);
        for t = 1:n.steps+1
            step_errors_bilin(t) = norm(comparison.q_BILINEAR(1:n.states, t, i) - comparison.q_nonlinear(:, t, i), 2);
            step_errors_edmd(t)  = norm(comparison.q_EDMD(1:n.states, t, i)    - comparison.q_nonlinear(:, t, i), 2);
            step_errors_lin(t)   = norm(comparison.q_LINEAR(:, t, i)            - comparison.q_nonlinear(:, t, i), 2);
        end
        err_bilin(i) = mean(step_errors_bilin);
        err_edmd(i)  = mean(step_errors_edmd);
        err_lin(i)   = mean(step_errors_lin);
    end
    
    comparison.errors.BILINEAR_2norm = mean(err_bilin);
    comparison.errors.EDMD_2norm     = mean(err_edmd);
    comparison.errors.LINEAR_2norm   = mean(err_lin);
    
    % 2. Testing errors on lifted states
    err_bilin_lifted = zeros(1, n.trajs_testing);
    err_edmd_lifted  = zeros(1, n.trajs_testing);
    
    for i = 1:n.trajs_testing
        step_errors_bilin = zeros(n.steps+1, 1);
        step_errors_edmd  = zeros(n.steps+1, 1);
        for t = 1:n.steps+1
            lifted_true          = f_lifting(comparison.q_nonlinear(:, t, i));  % lift column by column
            step_errors_bilin(t) = norm(comparison.q_BILINEAR(:, t, i) - lifted_true, 2);
            step_errors_edmd(t)  = norm(comparison.q_EDMD(:, t, i)     - lifted_true, 2);
        end
        err_bilin_lifted(i) = mean(step_errors_bilin);
        err_edmd_lifted(i)  = mean(step_errors_edmd);
    end
    
    comparison.errors.BILINEAR_2norm_lifted = mean(err_bilin_lifted);
    comparison.errors.EDMD_2norm_lifted     = mean(err_edmd_lifted);

    %4. Training erros
    % comparison.errors.BILINEAR_2norm_training = M_BILINEAR.training_error_2norm;
    comparison.errors.BILINEAR_fro_training = M_BILINEAR.training_error_fro;
    % comparison.errors.EDMD_2norm_training = M_EDMD.training_error_2norm;
    comparison.errors.EDMD_fro_training = M_EDMD.training_error_fro;
    % comparison.errors.LINEAR_2norm_training = M_LINEAR.training_error_2norm;
    comparison.errors.LINEAR_fro_training = M_LINEAR.training_error_fro;

    disp('Models compared');
end