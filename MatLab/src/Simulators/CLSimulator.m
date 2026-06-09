classdef CLSimulator < handle

    properties (Access=public)

        % System dynamics to simulate
        sys

        % FLController or SFController
        controller


        % Valid trajectory function
        discard_func

    end

    methods (Access=public)

        % Constructor
        function obj = CLSimulator(sys,controller,discard_func)

            obj.sys = sys;
            obj.controller = controller;

            if nargin < 3
                obj.discard_func = @(x,x0)(false);
            else
                obj.discard_func = @(x,x0)(discard_func(x,x0));
            end

        end

        function dataset = generateDataset(obj,nTrain,nTest,nSteps,initialStateFcn,q_ref,method)
        
            dataset.training = obj.generateTrajs(nTrain,nSteps,initialStateFcn,q_ref,method);        
            dataset.testing = obj.generateTrajs(nTest,nSteps,initialStateFcn,q_ref,method);
                
        end

        function trajs = generateTrajs(obj, nTrajs, nSteps, initialStateFcn, Ref, method)
                
            Xtot = zeros(obj.sys.nx, nSteps, nTrajs);
            Utot = zeros(obj.sys.nu, nSteps, nTrajs);
        
            for i = 1:nTrajs
                x0 = initialStateFcn();
                Xtot(:,1,i) = x0;
        
                for k = 1:nSteps-1
                    x_k = Xtot(:,k,i);
                    q_ref = Ref(:,k,i);
                    u_k = obj.controller.compute_u(x_k, q_ref);
                    Utot(:,k,i) = u_k;
        
                    switch method
                        case 'euler'
                            Xtot(:,k+1,i) = obj.sys.stepEuler(x_k, u_k);
                        case 'rk4'
                            Xtot(:,k+1,i) = obj.sys.stepRK4(  x_k, u_k);
                    end
        
                    if obj.discard_func(Xtot(:,k+1,i), x0)
                        break
                    end
                end
            end
        
            trajs.X = Xtot;
            trajs.U = Utot;
            trajs.Ref = Ref;

        end

    end

end