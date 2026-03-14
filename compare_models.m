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
            % q_lifted = f_lifting(q(:,j,i,2));
            q_edmdbilin(:,j+1,i) = M_BILINEAR.A*q_edmdbilin(:,j,i) + M_BILINEAR.B*u(:,j,i) + M_BILINEAR.N*q_edmdbilin(:,j,i)*u(:,j,i);
            % q(:,j+1,i,2) = M_BILINEAR.C*q_lifted_plus;
            
            % Linear EDMD model (lifting, computing evolution then going back to original states)
            % q_lifted = f_lifting(q(:,j,i,3));
            q_edmdlin(:,j+1,i) = M_EDMD.A*q_edmdlin(:,j,i) + M_EDMD.B*u(:,j,i);
            % q(:,j+1,i,3) = M_EDMD.C*q_lifted_plus;
            
            % Linear model
            q_lin(:,j+1,i) = M_LINEAR.A*q_lin(:,j,i) + M_LINEAR.B*u(:,j,i);
        end
    end
    comparison.q_nonlinear = q_nonlin(:,:,:);
    comparison.q_BILINEAR = q_edmdbilin(:,:,:);
    comparison.q_EDMD = q_edmdlin(:,:,:);
    comparison.q_LINEAR = q_lin(:,:,:);
    comparison.u = u;
    comparison.error_BILINEAR = norm(comparison.q_BILINEAR(1:n.states,:,:) - comparison.q_nonlinear, 'fro');
    comparison.error_EDMD = norm(comparison.q_EDMD(1:n.states,:,:) - comparison.q_nonlinear, 'fro');
    comparison.error_LINEAR = norm(comparison.q_LINEAR - comparison.q_nonlinear, 'fro');
    disp('Models compared');
end