classdef RobustController < handle

    properties (Access=public)

        FLModel
        K           % Controller (ss object)
        gamma       % Achieved H-inf norm
        method      % 'mixsyn' or 'datadriven'
        
        % Weights
        W1          % Performance weight
        W2          % Robustness weight  
        W3          % Input weight
        
        G_Set       % Multimodel set
        G_nominal
        relerr      % Errors wrt nominal model

        K_TF_DT     % Discrete transfer function (for LabVIEW)
        K_num
        K_den
        e_past
        v_past
        xK

    end
    
    methods (Access=public)
        
        function obj = RobustController(FLModel)

            obj.FLModel = FLModel;
        
        end
        
        function [u,v] = compute_u(obj, x, x_ref)

            if isempty(obj.xK)
                obj.xK = zeros(size(obj.K.A,1),1);
            end

            err = obj.FLModel.T(x) - obj.FLModel.T(x_ref);
            e = err(1);

            v = obj.K.C*obj.xK + obj.K.D*e;
            obj.xK = obj.K.A*obj.xK + obj.K.B*e;

            u = obj.FLModel.u(v, x);

        end
        
        function buildMultimodelSet(obj, data_robust, nfolds, Ts)
            
            ntrajs_training = size(data_robust.training.X, 3);
            traj_per_fold = floor(ntrajs_training / nfolds);
            
            G_set = [];
            
            for i = 1:nfolds
                y_fold = [];
                v_fold = [];
                for j = traj_per_fold*(i-1)+1 : traj_per_fold*i
                    for k = 1:size(data_robust.training.X, 2)
                        x_k = data_robust.training.X(:, k, j);
                        z_k = obj.FLModel.T(x_k);
                        y_fold(end+1, :) = z_k(1);
                    end
                    v_fold = [v_fold; squeeze(data_robust.training.V(1, :, j))'];
                end
                
                data_id = iddata(y_fold, v_fold, Ts);
                G_fold = ss(ssest(data_id, 2, 'Ts', Ts));

                if isempty(G_set)
                    G_set = G_fold;
                else
                    G_set = stack(1, G_set, G_fold);
                end

            end

            y_all = [];
            v_all = [];
            for j = 1:ntrajs_training
                y_all = [y_all; squeeze(data_robust.training.X(1, :, j))'];
                v_all = [v_all; squeeze(data_robust.training.V(1, :, j))'];
            end
            obj.G_nominal = ss(ssest(iddata(y_all, v_all, Ts), 2, 'Ts', Ts));            
            obj.G_Set = G_set;
            obj.relerr = obj.G_nominal\(obj.G_nominal-obj.G_Set);

            G_nominal_theory = c2d(ss(obj.FLModel.A_c, obj.FLModel.B_c, [1, 0], 0), Ts);

            figure;
            bodemag(obj.relerr,'b--', {0.1, pi/Ts});
            legend('Relerr', 'Location', 'southeast');
            title('Models error wrt nominal model');
            grid on;

            figure;
            bodemag(obj.G_Set,'b--', G_nominal_theory, 'r', {0.1, pi/Ts});
            legend('Identified models', 'Theoretical model', 'Location', 'southeast');
            title('Identified models wrt theoretical model');
            grid on;
            
        end
        
        function designMixsyn(obj, W1, W2, W3)
            obj.method = 'mixsyn';
            obj.W1 = W1;
            obj.W2 = W2;
            obj.W3 = W3;
            [obj.K, ~, obj.gamma] = mixsyn(obj.G_nominal, obj.W1, obj.W3, obj.W2);
            obj.computeTF();
            fprintf('Mixsyn done. Gamma = %.4f\n', obj.gamma);
        end
        
        function designDataDriven(obj, W1, W2, W3, Ts, ncontrolled_states, ninputs)
            obj.method = 'datadriven';
            obj.W1 = W1;
            obj.W2 = W2;
            obj.W3 = W3;
            
            % Frequency grid
            W2_frd = frd(obj.W2, logspace(-2, log10(pi/Ts), 50));
            omegas = unique([logspace(log10(0.3), log10(W2_frd.Frequency(end)), 200), ...
                             linspace(1, 10, 101), ...
                             linspace(100, 380, 101)]);
            
            % Augmented plant
            G_aug = augw(obj.G_nominal, obj.W1, obj.W3, obj.W2);
            G_aug = mktito(G_aug, 1, ninputs);
            
            % Initial stabilizing controller
            K_init = tf(0.001, 1, Ts);
            
            % Controller object
            K_ctrl = datadriven.Controller.SS(ncontrolled_states, 1, ninputs, Ts);
            K_ctrl.setinitial(K_init);
            
            % Synthesiser
            synth = Synthesiser(K_ctrl);
            synth.add_Hinf_objective(G_aug, omegas);
            synth.ensure_controller_stability(omegas);
            
            output = synth.synthesise();
            obj.K = output.Controller;
            obj.gamma = output.Objective;
            obj.computeTF();
            fprintf('Data-driven done.Gamma = %.4f\n', obj.gamma);
        end
        
        function validate(obj)
            nmodels = size(obj.G_Set, 3);
            fprintf('Validation across models \n');
            G_i = obj.G_nominal;
            S_i = feedback(1, G_i * obj.K);
            stable_S = isstable(S_i);
            fprintf('Nominal model : stable: %d \n', stable_S);
            for i = 1:nmodels
                G_i = obj.G_Set(:,:,i);
                S_i = feedback(1, G_i * obj.K);
                stable_S = isstable(S_i);
                fprintf('Fold %d : stable: %d \n', i, stable_S);
            end
        end
        
        function plotTFs(obj)
            S = feedback(1, obj.G_nominal * obj.K);
            T = feedback(obj.G_nominal * obj.K, 1);
            U = feedback(obj.K, obj.G_nominal);
            
            figure;
            bodemag(S, 1/obj.W1);
            legend('S(z)', '1/W1(z)', 'Location', 'southeast');
            title(['Sensitivity S - ' obj.method]);
            grid on;
            
            figure;
            bodemag(T, 1/obj.W2);
            legend('T(z)', '1/W2(z)', 'Location', 'southeast');
            title(['Complementary T - ' obj.method]);
            grid on;
            
            figure;
            bodemag(U, 1/obj.W3);
            legend('U(z)', '1/W3(z)', 'Location', 'southeast');
            title(['Input sensitivity U - ' obj.method]);
            grid on;

            figure;
            step(U);
            title(['Step Response - Control Signal (U) - ' obj.method]);
            grid on;
            
            figure;
            step(T);
            title(['Step response - Output (T) ' obj.method]);
            grid on;

        end
                
        function computeTF(obj)
            [b, a] = ss2tf(obj.K.A, obj.K.B, obj.K.C, obj.K.D);
            K_ct = tf(b, a);
            obj.K_TF_DT = c2d(K_ct, obj.K.Ts, 'zoh');
        
            obj.K_num = b(:)';
            obj.K_den = a(:)';
        
            nb = length(obj.K_num);
            na = length(obj.K_den);
        
            obj.e_past = zeros(1,nb);
            obj.v_past = zeros(1,na-1);
        
        end
        
    end
    
end