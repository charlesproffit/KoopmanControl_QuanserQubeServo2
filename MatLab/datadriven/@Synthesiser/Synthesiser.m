classdef Synthesiser < handle
    %Synthesiser Frequency-domain data-driven controller synthesiser
    % 
    %   Synthesiser(controller) gives a optimiser object for.
    % 
    %   Synthesiser(..., "verbose", verbosity) sets the verbosity. 
    %       (Default: true)
    % 
    %   Synthesiser(..., "solver", solver) set the synthesiser to use the
    %   `solver` for optimisation. 
    %       (Default: "mosek")
    % 
    %   ------------------------------------------------------------------------
    %   Copyright 2025 Vaibhav Gupta, DDMAC, EPFL (MIT License)
    %            

    properties (SetAccess=protected)
        H2_constraints = {}     % H2 Constraints
        Hinf_constraints = {}   % Hinf Constraints

        H2_objectives = {}      % H2 Objectives
        Hinf_objectives = {}    % Hinf Objectives

        controller_stability = {} % Ensure constroller stability

        controller  % Controller
    end

    properties (Access=protected)
        logger      % Logger
        sdp_opts    % Optimisation options
    end

    %% Constructor Method
    methods
        function obj = Synthesiser(controller, opts)
            arguments
                controller (1, 1) datadriven.Controller.Controller
                opts.verbose (1, 1) logical = true
                opts.solver  (1, 1) string = "mosek"
            end
            
            obj.controller = controller;

            % Initialize optimisation logger
            fid = 1;
            if ~opts.verbose
                fid = 0;
            end
            obj.logger = utils.Logger(fid, ...
                "Slack"     , "% 10.4f", ...
                "Objective" , "% 10.4f", ...
                "Decrease"  , "% 10.4f", ...
                "Total Time", "% 10.3f"  , ...
                "Solve Time", "% 10.3f"  , ...
                "Status"    , "%-40s"   );
            
            % Optimisation Settings
            obj.sdp_opts = sdpsettings('solver', opts.solver, 'verbose', 0);
        end
    end

    %% Add objectives and constraints
    methods
        function ensure_controller_stability(obj, omegas)
            arguments
                obj
                omegas (:, 1) double
            end
            obj.controller_stability = { obj.controller, omegas };
        end

        function add_H2_objective(obj, plants, omegas)
            arguments
                obj
                plants lti
                omegas (:, 1) double
            end
            n_plants = size(plants, 3);
            
            for i_plant = 1:n_plants
                obj.H2_objectives = [obj.H2_objectives;
                    {{ plants(:, :, i_plant), obj.controller, omegas }};
                    ];
            end
        end

        function add_Hinf_objective(obj, plants, omegas)
            arguments
                obj
                plants lti
                omegas (:, 1) double
            end
            n_plants = size(plants, 3);
            
            for i_plant = 1:n_plants
                obj.Hinf_objectives = [obj.Hinf_objectives;
                    {{ plants(:, :, i_plant), obj.controller, omegas }};
                    ];
            end
        end
    
        function add_H2_constraint(obj, plants, omegas)
            arguments
                obj
                plants lti
                omegas (:, 1) double
            end
            n_plants = size(plants, 3);
            
            for i_plant = 1:n_plants
                obj.H2_constraints = [obj.H2_constraints;
                    {{ plants(:, :, i_plant), obj.controller, omegas }};
                    ];
            end
        end

        function add_Hinf_constraint(obj, plants, omegas)
            arguments
                obj
                plants lti
                omegas (:, 1) double
            end
            n_plants = size(plants, 3);
            
            for i_plant = 1:n_plants
                obj.Hinf_constraints = [obj.Hinf_constraints;
                    {{ plants(:, :, i_plant), obj.controller, omegas }};
                    ];
            end
        end
    end
    
end