classdef OLSimulator < handle

    properties (Access=public)

        % System dynamics to simulate
        sys

        % Valid trajectory function
        discard_func

    end

    methods (Access=public)

        % Constructor
        function obj = OLSimulator(sys,discard_func)

            obj.sys = sys;
            if nargin < 2
                obj.discard_func = @(x,x0)(false);
            else
                obj.discard_func = @(x,x0)(discard_func(x,x0));
            end

        end

        function X = simulate_onestep(obj,Xreal,U,method)
        
            nsteps = size(U,2);
            ntrajs = size(U,3);

            if ismethod(obj.sys, 'liftState')
                X = zeros(obj.sys.nz, nsteps, ntrajs);
            else
                X = zeros(obj.sys.nx, nsteps, ntrajs);
            end
        
            for i=1:ntrajs
                if ismethod(obj.sys, 'liftState')
                    x0 = obj.sys.liftState(Xreal(:,1,i));
                else
                    x0 = Xreal(:,1,i);
                end
                 X(:,1,i) = x0;
    
                for k=1:nsteps-1
    
                    if ismethod(obj.sys, 'liftState')
                        x_k = obj.sys.liftState(Xreal(:,k,i));
                    else
                        x_k = Xreal(:,k,i);
                    end
            
                    switch method
            
                        case 'euler'
                            X(:,k+1,i) = obj.sys.stepEuler(x_k,U(:,k,i));
            
                        case 'rk4'
                            X(:,k+1,i) = obj.sys.stepRK4(x_k,U(:,k,i));
            
                    end
            
                end

            end
            
        end

        function X = simulate_multistep(obj,X0,U,method)
        
            nsteps = size(U,2);
            ntrajs = size(U,3);

            if ismethod(obj.sys, 'liftState')
                X = zeros(obj.sys.nz, nsteps, ntrajs);
            else
                X = zeros(obj.sys.nx, nsteps, ntrajs);
            end

            for i=1:ntrajs
                if ismethod(obj.sys, 'liftState')
                    x0 = obj.sys.liftState(X0(:,i));
                else
                    x0 = X0(:,i);
                end
                 X(:,1,i) = x0;
            
                for k=1:nsteps-1
            
                    switch method
            
                        case 'euler'
                            X(:,k+1,i) = obj.sys.stepEuler(X(:,k,i),U(:,k,i));
            
                        case 'rk4'
                            X(:,k+1,i) = obj.sys.stepRK4(X(:,k,i),U(:,k,i));
            
                    end
            
                end

            end
        
        end

        function [X,valid] = simulate_with_check(obj,x0,U,method)
        
            valid = true;

            nsteps = size(U,2);
        
            X = zeros(obj.sys.nx,nsteps);
            X(:,1) = x0;
        
            for k=1:nsteps-1
        
                switch method
        
                    case 'euler'
                        X(:,k+1) = obj.sys.stepEuler(X(:,k),U(:,k));
        
                    case 'rk4'
                        X(:,k+1) = obj.sys.stepRK4(X(:,k),U(:,k));
        
                end

                if obj.discard_func(X(:,k+1), x0)

                    valid = false;
                    X = X(:,1:k+1);
                    return

                end
        
            end
        
        end

        function dataset = generateDataset(obj,nTrain,nTest,nSteps,initialStateFcn,inputGeneratorFcn,method)
        
            dataset.training = obj.generateTrajs(nTrain,nSteps,initialStateFcn,inputGeneratorFcn, method);        
            dataset.testing = obj.generateTrajs(nTest,nSteps,initialStateFcn,inputGeneratorFcn,method);
                
        end

        function trajs = generateTrajs(obj,nTrajs,nSteps,initialStateFcn,inputGeneratorFcn,method)
    
            Xtot = zeros(obj.sys.nx,nSteps,nTrajs);
            Utot = zeros(obj.sys.nu,nSteps,nTrajs);
        
            for i = 1:nTrajs
    
                valid = false;
    
                while ~valid
    
                    x0 = initialStateFcn();
    
                    U = inputGeneratorFcn();
    
                    [X,valid] = obj.simulate_with_check(x0,U,method);
    
                    if valid
    
                        Xtot(:,:,i) = X;
                        Utot(:,:,i) = U;
    
                    end
    
                end
    
            end
    
            trajs.X = Xtot;
            trajs.U = Utot;

        end

    end

end