classdef EDMDModel < handle

    properties (Access=public)
    
        % Identification settings
        mode           % 'LINEAR','BILINEAR'
        method         % 'LS','RIDGE'
    
        % Timestep
        Ts

        % Dictionary
        f_lifting
    
        % Dimensions
        nx
        nz
        nu
        ncv
        nsteps
        ntrajs

        % Identified model
        A_CT
        B_CT
        A_DT
        B_DT
        C
    
        % Statistics
        lambda
        condition_number
        svds
    
        training_error_fro
        training_error_fro_lifted
    
    end

    methods (Access=public)

        function obj = EDMDModel(mode,method,f_lifting,nx,nu,nz,Ts,ncrossval)
        
            obj.mode = mode;
            obj.method = method;
        
            obj.f_lifting = @(q)f_lifting(q);
        
            obj.nx = nx;
            obj.nu = nu;
            obj.nz = nz;
            obj.ncv = ncrossval;

            obj.Ts = Ts;
        
        end

        function EDMDIdentification(obj,dataset)
            
            [Z,Z_plus,Z_dot] = obj.buildMatrices(dataset);

            [obj.A_CT, obj.B_CT] = obj.identifyModel(Z,Z_dot);
            [obj.A_DT, obj.B_DT] = obj.identifyModel(Z,Z_plus);

        end

        function [Z,Z_plus,Z_dot] = buildMatrices(obj,dataset)

            obj.nsteps = size(dataset.training.X,2);
            obj.ntrajs = size(dataset.training.X,3);

            % Step 1 : Compute lifted states
            u = dataset.training.U;
            q = dataset.training.X;
            phi = zeros(obj.nz, obj.nsteps, obj.ntrajs);
            for i = 1:obj.ntrajs
                for j = 1:obj.nsteps
                    phi(:,j,i) = obj.f_lifting(q(:,j,i));
                end
            end
        
            % Step 2 : Stack data
            U = reshape(u, obj.nu, []);
            Q = reshape(q, obj.nx, []);
            Phi  = reshape(phi, obj.nz, []);
            PhiU = Phi .* repmat(U, obj.nz, 1);

            obj.C = Q * pinv(Phi); %% C miniminizes ||Q-C*Phi||=> C goes from lifted states to original states

            normal_idx = [];
            plus_idx  = [];
            for i = 1:obj.ntrajs
                base = (i-1)*obj.nsteps;
                normal_idx = [normal_idx, base + (1:obj.nsteps-1)];
                plus_idx  = [plus_idx,  base + (2:obj.nsteps)];
            end

            U_normal = U(:,normal_idx);
            U_plus = U(:,plus_idx);
            Phi_normal = Phi(:,normal_idx);
            Phi_plus = Phi(:,plus_idx);
            PhiU_normal = PhiU(:,normal_idx);
            PhiU_plus = PhiU(:,plus_idx);

            if obj.mode == "BILINEAR"
                Z = [Phi_normal ; PhiU_normal];
                Z_plus = [Phi_plus ; PhiU_plus];
            else
                Z = [Phi_normal ; U_normal];
                Z_plus = [Phi_plus ; U_plus];
            end
            
            obj.svds = svd(Z_plus);
            obj.condition_number = cond(Z_plus);
            Z_dot = (Z_plus - Z) / obj.Ts;

        end

        function [A,B] = identifyModel(obj,Z,Z_plus)

            if obj.method == "RIDGE"
                obj.lambda = obj.selectLambda(Z, Z_plus);
                AB = Z_plus * Z' * (Z*Z' + obj.lambda*eye(size(Z,1)))^(-1);
            else
                obj.lambda = 0;
                AB = Z_plus * pinv(Z);
            end
            A = AB(1:obj.nz, 1:obj.nz);
            B = AB(1:obj.nz, obj.nz+1:end);

            % Training error
            diff = Z_plus - AB*Z;
            diff = diff(1:obj.nz, :);
            lifted_diff = reshape(diff,obj.nz,obj.nsteps-1, obj.ntrajs);
            state_diff = obj.C * diff;
            state_diff = reshape(state_diff, obj.nx, obj.nsteps-1, obj.ntrajs);
            
            traj_fro= zeros(obj.ntrajs,1);
            traj_fro_lifted = zeros(obj.ntrajs,1);
            for i = 1:obj.ntrajs
                traj_fro(i) = norm(state_diff(:,:,i), 'fro');
                traj_fro_lifted(i) = norm(lifted_diff(:,:,i), 'fro');
            end
            obj.training_error_fro = mean(traj_fro);
            obj.training_error_fro_lifted = mean(traj_fro_lifted);

        end

        function lambda = selectLambda(obj, Z, Z_target)
            lambdas = [0, logspace(-4, 4, 50)];
            n_lambdas = length(lambdas);

            traj_fold_assignments = mod(0:obj.ntrajs-1, obj.ncv) + 1;
            fold_assignments = repelem(traj_fold_assignments, obj.nsteps-1);

            val_errors = zeros(obj.ncv, n_lambdas);

            for fold = 1:obj.ncv
                val_mask = (fold_assignments == fold);
                train_mask = ~val_mask;

                Z_train = Z(:, train_mask);
                Z_val = Z(:, val_mask);
                target_train = Z_target(:, train_mask);
                target_val = Z_target(:, val_mask);

                for k = 1:n_lambdas
                    if lambdas(k) == 0
                        AB = target_train * pinv(Z_train);
                    else
                        AB = target_train * Z_train' * (Z_train*Z_train' + lambdas(k)*eye(size(Z_train,1)))^(-1);
                    end
                    val_errors(fold, k) = norm(target_val - AB*Z_val, 'fro');
                end
            end

            mean_val_errors = mean(val_errors, 1);
            [~, best_k] = min(mean_val_errors);
            lambda = lambdas(best_k);
            fprintf('[%s/%s] Best lambda = %.6f\n', obj.mode, obj.method, lambda);

            % figure;
            % lambda_plot    = lambdas;
            % lambda_plot(1) = 1e-5;
            % semilogx(lambda_plot, mean_val_errors, 'b-o', 'LineWidth', 2, 'MarkerSize', 4);
            % xline(max(lambda, 1e-5), 'r--', 'LineWidth', 2, ...
            %       'DisplayName', sprintf('\\lambda^* = %.2e', lambda));
            % xlabel('\lambda');
            % ylabel('Mean validation error (Frobenius norm)');
            % title(sprintf('Cross-validation curve - %s/%s', obj.mode, obj.method));
            % grid on;
        end

        % Continuous dynamics
        function xdot = continuousDynamics(obj,x,u)
            if obj.mode == "BILINEAR"
                xdot = obj.A_CT*x + obj.B_CT*(x.*u);
            else
                xdot = obj.A_CT*x + obj.B_CT*u;
            end
        end

        % Euler step
        function xplus = stepEuler(obj,x,u)

            xplus = x + obj.Ts*obj.continuousDynamics(x,u);

        end

        % RK4 step
        function xplus = stepRK4(obj,x,u)

            k1 = obj.continuousDynamics(x,u);

            k2 = obj.continuousDynamics( x + obj.Ts/2*k1 , u);

            k3 = obj.continuousDynamics( x + obj.Ts/2*k2 , u);

            k4 = obj.continuousDynamics( x + obj.Ts*k3 , u);

            xplus = x + obj.Ts/6*( k1 + 2*k2 + 2*k3 + k4 );

        end

        function z = liftState(obj, x)
            z = obj.f_lifting(x);
        end

        function x = projectState(obj, z)
            x = obj.C*z;
        end

    end

end