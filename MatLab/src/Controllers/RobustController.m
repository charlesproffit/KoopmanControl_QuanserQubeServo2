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
        G_nominal_theory    % Nominal Brunovsky matrices model
        relerr_theory

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
            
            function buildMultimodelSet(obj, data_robust, nfolds, Ts)
                
                ntrajs_training = size(data_robust.training.X, 3);
                traj_per_fold = floor(ntrajs_training / nfolds);
                
                G_set = [];
    
                freqs = logspace(log10(0.1),log10(0.999*pi/Ts),400);

                for i = 1:nfolds
                    y = [];
                    u = [];
                    for j = traj_per_fold*(i-1)+1 : traj_per_fold*i
                        y = [y; data_robust.training.X(1, :, j)'];
                        u = [u; data_robust.training.V(:, :, j)'];
                    end
                    if isempty(G_set)
                        G_set = spa(iddata(y, u, Ts), [], freqs);
                    else
                        G_set = stack(1, G_set, spa(iddata(y, u, Ts), [], freqs));
                    end
                end
                obj.G_Set = G_set;
    
                obj.G_nominal_theory = c2d(ss(obj.FLModel.A_c, obj.FLModel.B_c, [1, 0], 0), Ts);
                obj.relerr_theory = obj.G_nominal_theory\(obj.G_nominal_theory-obj.G_Set);
                
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
            G = ss(obj.G_nominal);
            S = feedback(1, G * obj.K);
            T = feedback(G * obj.K, 1);
            
            [gm, pm] = margin(G * obj.K);
            stable = isstable(S);
            
            fprintf('Nominal model validation:\n');
            fprintf('Stable: %d | GM=%.2f dB | PM=%.2f deg\n', ...
                stable, 20*log10(gm), pm);
            fprintf('Peak S: %.2f dB | Peak T: %.2f dB\n', ...
                20*log10(norm(S,'inf')), 20*log10(norm(T,'inf')));
        end
    
        function plotTFs(obj)
            G = ss(obj.G_nominal);
            S = feedback(1, G * obj.K);
            T = feedback(G * obj.K, 1);
            U = feedback(obj.K, G);
            
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