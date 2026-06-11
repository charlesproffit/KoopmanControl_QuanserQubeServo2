classdef SS < datadriven.Controller.Controller
    %SS State-space Controller object
    % 
    %   SS(order, nMeasurements, nActuators) gives a continuous-time
    %   State-space Controller object.
    % 
    %   SS(order, nMeasurements, nActuators, sampleTime) gives a 
    %   State-space Controller object with given sample time.
    % 
    %   SS(order, nMeasurements, nActuators, sampleTime, factorisation) 
    %   gives a State-space Controller object with selected factorisation
    %   ("Left" or "Right").
    % 
    %   ------------------------------------------------------------------------
    %   Copyright 2025 Vaibhav Gupta, DDMAC, EPFL (MIT License)
    %
    
    %% Constructor Method
    methods
        function obj = SS(order, nMeasurements, nActuators, sampleTime, factorisation)
            arguments (Input)
                order           (1, 1) double
                nMeasurements   (1, 1) double
                nActuators      (1, 1) double
                sampleTime      (1, 1) double = 0
                factorisation   (1, 1) datadriven.Controller.Factorisation = "Left"
            end
            
            obj.nx = order;
            obj.ny = nMeasurements;
            obj.nu = nActuators;
            obj.Ts = sampleTime;
            obj.factorisation = factorisation;
            
            obj.initialise();
        end
    end

    %% Methods
    methods
        function setinitial(obj, sys)
            %SETINITIAL Sets the initial controller for the optimisation
            % 
            %   This function is designed to set the initial stabilising
            %   controller.
            % 
            arguments (Input)
                obj 
                sys lti
            end
            
            [nu_sys, ny_sys] = size(sys);
            if (nu_sys ~= obj.nu) || (ny_sys ~= obj.ny)
                error("DataDriven:SS:MismatchIO", ...
                    "Number of input/output is not correct!");
            end
            
            nx_sys = order(sys);
            if (nx_sys > obj.nx)
                error("DataDriven:SS:MismatchOrder", ...
                    "Cannot use higher order controller as the initial!");
            end
            
            obj.initialise(ss(sys));
        end

        function normalise(obj)
            %NORMALISE Normalise the controller for better numerical stability
            %   
            arguments (Input)
                obj 
            end
            obj.setinitial(balreal(ss(obj)));
        end
    end

    %% Overloaded methods 
    methods        
        function [X, Y] = freqresp(obj, omegas) 
            %FREQRESP Evaluates frequency response over a grid of frequencies
            % 
            %   This function evaluates the frequency response of the SS object
            %   at given frequency points, 'omegas' rad/s.
            % 
            %   Outputs are returned as frequency response of X and frequency
            %   response of Y which are defined as,
            %       * Y\X : for 'Left' factorisation
            %       * X/Y : for 'Right' factorisation
            % 
            %   FREQRESP is a more complex version of EVALFR meant for
            %   evaluation of the response over a grid of frequencies. Use
            %   EVALFR to compute the frequency response at a single point.
            % 
            arguments (Input)
                obj 
                omegas (1, 1, :) double
            end
            
            % ----- Internally in Yalmip is the loop too! -----
            % if obj.Ts == 0
            %     z = 1j * omegas;
            % elseif obj.Ts == -1
            %     z = exp(1j * omegas);
            % else
            %     z = exp(1j * omegas * obj.Ts);
            % end
            % 
            % switch (obj.factorisation)
            %     case {"Right"}
            %         XY = pagemtimes(obj.C, pagemldivide(eye(obj.nx).*z - obj.A, obj.B)) + obj.D;
            %         Y = XY(1:obj.ny, :, :);
            %         X = XY(obj.ny+1:end, :, :);
            %     case {"Left"}
            %         YX = pagemtimes(pagemrdivide(obj.C, eye(obj.nx).*z - obj.A), obj.B) + obj.D;
            %         Y = YX(:, 1:obj.nu, :);
            %         X = YX(:, obj.nu+1:end, :);
            %     otherwise
            %         error("datadriven:Controller:InvalidFactorisation", ...
            %             "Factorisation type is invalid!")
            % end

            % TODO: Write a vectorised version of FREQRESP.
            nOmegas = length(omegas);

            if nOmegas == 1
                [X, Y] = evalfr(obj, omegas);
            else
                X = cell(1, nOmegas);
                Y = cell(1, nOmegas);
                for iOmegas = 1:nOmegas
                    [X{iOmegas}, Y{iOmegas}] = obj.evalfr(omegas(iOmegas));
                end
                % Reaches max recursion if `nOmegas > ~500`
                % X = cat(3, X{:});
                % Y = cat(3, Y{:});

                X = reshape([X{:}], [ size(X{1}) , nOmegas ]);
                Y = reshape([Y{:}], [ size(Y{1}) , nOmegas ]);
            end
        end

        function [X, Y] = evalfr(obj, omega)
            %EVALFR Evaluates frequency response at a single frequency
            % 
            %   This function evaluates the frequency response of the SS object
            %   at a single frequency point, 'omega' rad/s.
            % 
            %   Outputs are returned as frequency response of X and frequency
            %   response of Y which are defined as,
            %       * Y\X : for 'Left' factorisation
            %       * X/Y : for 'Right' factorisation
            % 
            %   EVALFR is a simplified version of FREQRESP meant for quick
            %   evaluation of the response at a single point. Use FREQRESP to
            %   compute the frequency response over a grid of frequencies.
            % 
            arguments (Input)
                obj
                omega (1, 1) double
            end

            if obj.Ts == 0
                z = 1j * omega;
            elseif obj.Ts == -1
                z = exp(1j * omega);
            else
                z = exp(1j * omega * obj.Ts);
            end
            
            switch (obj.factorisation)
                case {"Right"}
                    XY = obj.C / (eye(obj.nx)*z - obj.A) * obj.B + obj.D;
                    Y = XY(1:obj.ny, :);
                    X = XY(obj.ny+1:end, :);
                case {"Left"}
                    YX = obj.C / (eye(obj.nx)*z - obj.A) * obj.B + obj.D;
                    Y = YX(:, 1:obj.nu);
                    X = YX(:, obj.nu+1:end);
                otherwise
                    error("datadriven:Controller:InvalidFactorisation", ...
                        "Factorisation type is invalid!")
            end
        end
       
        function K = ss(obj)
            %SS Convert from State-space Controller object to MATLAB SS object
            % 
            arguments (Input)
                obj 
            end

            switch (obj.factorisation)
                case {"Right"}
                    Cy = double(obj.C(1:obj.ny, :)); 
                    Cx = double(obj.C(obj.ny+1:end, :));
                    Dy = double(obj.D(1:obj.ny, :)); 
                    Dx = double(obj.D(obj.ny+1:end, :));
    
                    K = ss(...
                        obj.A - (obj.B/Dy)*Cy, ...
                        obj.B/Dy, ...
                        Cx - (Dx/Dy)*Cy, ...
                        Dx/Dy, ...
                        obj.Ts);
                case {"Left"}
                    By = double(obj.B(:, 1:obj.nu)); 
                    Bx = double(obj.B(:, obj.nu+1:end));
                    Dy = double(obj.D(:, 1:obj.nu)); 
                    Dx = double(obj.D(:, obj.nu+1:end));
    
                    K = ss(...
                        obj.A - By*(Dy\obj.C), ...
                        Bx - By*(Dy\Dx), ...
                        Dy\obj.C, ...
                        Dy\Dx, ...
                        obj.Ts);
                otherwise
                    error("datadriven:Controller:InvalidFactorisation", ...
                        "Factorisation type is invalid!")
            end
        end
    end
    
    %% Internal Methods
    methods (Access=protected)
        function initialise(obj, initialSystem)
            %INITIALISE Configures the internal and optimisation variables
            % 
            %   This function configures the object either with random 0 system
            %   or the given 'initialSystem'.
            % 
            arguments (Input)
                obj
                initialSystem lti = ss([], [], [], zeros(obj.nu, obj.ny), obj.Ts)
            end

            K = obj.augmentsystem(initialSystem, obj.nx, obj.factorisation);

            switch (obj.factorisation)
                case {"Right"}
                    [A_tmp, B_tmp, C_tmp, D_tmp, ~] = ssdata(rncf(balreal(K)));
                case {"Left"}
                    [A_tmp, B_tmp, C_tmp, D_tmp, ~] = ssdata(lncf(balreal(K)));
                otherwise
                    error("datadriven:Controller:InvalidFactorisation", ...
                        "Factorisation type is invalid!")
            end
            
            % ===== Assign SDPVARs =====
            obj.A = A_tmp;
            switch (obj.factorisation)
                case {"Right"}
                    obj.B = B_tmp;
                    obj.C = sdpvar(size(C_tmp, 1), size(C_tmp, 2), 'full');
                    if ~isempty(obj.C), assign(obj.C, C_tmp); end
                case {"Left"}
                    obj.B = sdpvar(size(B_tmp, 1), size(B_tmp, 2), 'full');
                    if ~isempty(obj.B), assign(obj.B, B_tmp); end
                    obj.C = C_tmp;
                otherwise
                    error("datadriven:Controller:InvalidFactorisation", ...
                        "Factorisation type is invalid!")
            end
            obj.D = sdpvar(size(D_tmp, 1), size(D_tmp, 2), 'full');
            assign(obj.D, D_tmp);
        end
    end

    methods (Static, Access=protected)
        function augmentedSystem = augmentsystem(originalSystem, desiredOrder, factorisation, options)
            %AUGMENTSYSTEM Augment system to desired order
            % 
            %   This function augments the 'originalSystem' by adding extra
            %   states to it in order to attain the desired order, while
            %   ensuring that the original system's dynamics remain unchanged.
            % 
            arguments (Input)
                originalSystem lti
                desiredOrder (1, 1) int32
                factorisation (1, 1) datadriven.Controller.Factorisation
                options.stabilityDistance (1, 1) double = 0.01
            end
            
            if (order(originalSystem) > desiredOrder)
                error("DataDriven:SS:augmentsystem:OrderTooHigh", ...
                    "Cannot use higher order controller as the initial!");
            end
            
            % ===== Zero Gain System =====
            nx_zero = desiredOrder - order(originalSystem);
            
            % Generate random A matrix with stable poles
            %   w/ 0.5 probability of complex and real roots
            n_complex = floor(sum(rand(nx_zero, 1) < 0.5) / 2);
            n_real = nx_zero - 2*n_complex;

            if originalSystem.Ts == 0      % Continuous time
                complex_eig = ...
                    - exp(randn(n_complex, 1)) - 1 - options.stabilityDistance ...
                    + 3j * exp(randn(n_complex, 1));
                real_eig = - exp(randn(n_real, 1)) - 1;
            else                % Discrete time
                complex_eig = (1 - options.stabilityDistance) * rand(n_complex, 1) .* exp(1j * pi * rand(n_complex, 1));
                real_eig = 2 * rand(n_real, 1) - 1;
            end
            tmp = [complex_eig, conj(complex_eig)].';
            A_eig_complex = diag([tmp(:); real_eig]);

            % Convert to real block diagonal form
            [~, A_eig] = cdf2rdf(eye(desiredOrder), A_eig_complex);

            % Random orthoganal basis multiplication
            T = orth(randn(nx_zero));
            A_tmp = T \ A_eig * T;
            
            [nOutputs, nInputs] = size(originalSystem);
            switch (factorisation)
                case {"Right"}
                    B_tmp = randn(nx_zero, nInputs);
                    C_tmp = zeros(nOutputs, nx_zero);
                case {"Left"}
                    B_tmp = zeros(nx_zero, nInputs);
                    C_tmp = randn(nOutputs, nx_zero);
                otherwise
                    B_tmp = zeros(nx_zero, nInputs);
                    C_tmp = zeros(nOutputs, nx_zero);
            end
            
            D_tmp = zeros(nOutputs, nInputs);
            
            % ===== Augmented System =====
            augmentedSystem = ss(...
                blkdiag(originalSystem.A, A_tmp), ...
                [originalSystem.B; B_tmp], ...
                [originalSystem.C, C_tmp], ...
                originalSystem.D + D_tmp, ...
                originalSystem.Ts);
        end
    end

    %% Properties
    properties (SetAccess=protected)
        % Order of the controller
        nx (1, 1) double
    end
    
    properties (Access=protected)
        A
        B
        C
        D
    end
    
end