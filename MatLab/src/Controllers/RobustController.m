classdef RobustController < handle

    properties (Access=public)

        FLModel
        K           % Controller (ss object)
        K_TF_DT     % Discrete transfer function (for LabVIEW)
        gamma       % Achieved H-inf norm
        method      % 'mixsyn' or 'datadriven'
        
        % Weights
        W1          % Performance weight
        W2          % Robustness weight  
        W3          % Input weight
        
        G_Set       % Multimodel set
        G_nominal
        relerr      % Errors wrt nominal model
    end
    
    methods (Access=public)
        
        function obj = RobustController(FLModel)

            obj.FLModel = FLModel;
        
        end
        
        function u = compute_u(obj, x, x_ref)

            err = obj.FLModel.T(x) - obj.FLModel.T(x_ref);
            v = -obj.K * err;
            u = obj.FLModel.u(v, x);

        end
        
        function buildMultimodelSet(obj, data_robust, nfolds, func_linear, ncontrolled_states, ninputs, Ts, ncross_val_groups)
            
            ntrajs_training = size(data_robust.training.X, 3);
            traj_per_fold = floor(ntrajs_training / nfolds);
            
            G_set = [];
            A_sum = zeros(size(obj.FLModel.A_c));
            B_sum = zeros(size(obj.FLModel.B_c));
            
            for i = 1:nfolds
                Model_ROBUST = EDMDModel("LINEAR", "LS", func_linear, ncontrolled_states, ninputs, ncontrolled_states, Ts, ncross_val_groups);
                
                dataset_fold.training.X = data_robust.training.X([1,3], :, traj_per_fold*(i-1)+1 : traj_per_fold*i);
                dataset_fold.training.U = data_robust.training.U(:, :, traj_per_fold*(i-1)+1 : traj_per_fold*i);
                
                Model_ROBUST.EDMDIdentification(dataset_fold);
                
                A_sum = A_sum + Model_ROBUST.A_CT;
                B_sum = B_sum + Model_ROBUST.B_CT;
                
                if isempty(G_set)
                    G_set = c2d(ss(Model_ROBUST.A_CT, Model_ROBUST.B_CT, [1, 0], 0), Ts);
                else
                    G_set = stack(1, G_set, c2d(ss(Model_ROBUST.A_CT, Model_ROBUST.B_CT, [1, 0], 0), Ts));
                end
            end
            
            obj.G_Set = G_set;
            obj.G_nominal = c2d(ss(A_sum / nfolds, B_sum / nfolds, [1, 0], 0), Ts);
            % obj.G_nominal = c2d(ss(obj.FLModel.A_c, obj.FLModel.B_c, [1, 0], 0), Ts);

            obj.relerr = obj.G_nominal\(obj.G_nominal-obj.G_Set);

            figure;
            bodemag(obj.relerr,'b--', {0.1, pi/Ts});
            legend('Relerr', 'Location', 'southeast');
            title('Multimodel set error wrt nominal model');
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
            step(T);
            title(['Step response - ' obj.method]);
            grid on;

        end
                
        function computeTF(obj)
            [b, a] = ss2tf(obj.K.A, obj.K.B, obj.K.C, obj.K.D);
            K_ct = tf(b, a);
            obj.K_TF_DT = c2d(K_ct, obj.K.Ts, 'zoh');
        end
        
    end
    
end