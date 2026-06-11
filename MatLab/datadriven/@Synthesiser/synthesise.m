function outputs = synthesise(obj, opts)
    arguments
        obj 
        opts.maxIter (1, 1) double = 20
    end

    constrain_satisfied = isempty(obj.Hinf_constraints) && isempty(obj.H2_constraints);
    
    obj.logger.init();
    prev_obj = Inf;
    for iter = 1:opts.maxIter
        obj.logger.log_iter(iter);
        t_begin = tic();
        
        LMIs = [];
        % Ensure controller stability
        if ~isempty(obj.controller_stability)
            LMIs = datadriven.Constraint.controller_stability( obj.controller_stability{:} );
        end

        % Create constraints
        slack = 0;

        for i = 1:length(obj.Hinf_constraints)
            [local_LMIs, local_slack] = datadriven.Constraint.Hinf( obj.Hinf_constraints{i}{:}, "strict", constrain_satisfied);
            slack = slack + local_slack;
            LMIs = [LMIs; local_LMIs];
        end

        for i = 1:length(obj.H2_constraints)
            [local_LMIs, local_slack] = datadriven.Constraint.H2( obj.H2_constraints{i}{:}, "strict", constrain_satisfied);
            slack = slack + local_slack;
            LMIs = [LMIs; local_LMIs];
        end
        
        % Create objectives
        if constrain_satisfied
            objective = sdpvar();

            for i = 1:length(obj.Hinf_objectives)
                [local_LMIs, local_objective] = datadriven.Objective.Hinf( obj.Hinf_objectives{i}{:} );
                LMIs = [LMIs; local_LMIs; local_objective <= objective];
            end

            for i = 1:length(obj.H2_objectives)
                [local_LMIs, local_objective] = datadriven.Objective.H2( obj.H2_objectives{i}{:} );
                LMIs = [LMIs; local_LMIs; local_objective <= objective];
            end
        end

        % Controller optimisation
        if constrain_satisfied
            JOB = optimize(LMIs, objective, obj.sdp_opts);
        else
            JOB = optimize(LMIs, slack, obj.sdp_opts);
        end
        
        if JOB.problem ~= 0
            % Quit on optimisation failure
            obj.logger.log_err(strcat(...
                "Optimization Failed! => ", ...
                JOB.info));
            break;
        end
        
        % Normalise and extract controller
        obj.controller.normalise();
        Kc = ss(obj.controller);

        total_time = toc(t_begin);

        if constrain_satisfied
            curr_obj = sqrt(double(objective));
            % Log the interation outputs
            obj.logger.log(...
                double(slack), ...
                curr_obj, ...
                curr_obj - prev_obj, ...
                total_time, ...
                JOB.solvertime, ...
                JOB.info);
            
            % Check for convergence
            err = abs(curr_obj - prev_obj);
            if (err < 1e-4) || (err/prev_obj < 1e-4)
                break;
            end
            prev_obj = curr_obj;
        else
            % Log the interation outputs
            obj.logger.log(...
                double(slack), ...
                [], ...
                [], ...
                total_time, ...
                JOB.solvertime, ...
                JOB.info);
            
            % Check for constrain satisfaction
            err = abs(double(slack));
            if (err < 1e-4)
                constrain_satisfied = true;
            end
        end
    end
    obj.logger.cleanup();
    
    % Output Structure
    outputs = struct( ...
        "Objective", curr_obj, ...
        "Controller", Kc ...
        );
end
