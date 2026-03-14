function [data,n] = collect_data(f_discrete, umax, n, mode)
    if strcmp(mode, 'EDMD')
        % 1. TRAINING DATA
        n.trajs_training = fix(n.train_test_ratio*n.trajs);
        angle_region = 2*pi/n.regions;
        q = zeros(n.states, n.steps+1, n.trajs_training);
        u = 2*umax*rand(n.inputs, n.steps, n.trajs_training) - umax; % [-umax, umax]
    
        % Reminder : States q :
        %   q(1) => Angle of rotary arm
        %   q(2) => Angle of pendulum arm
        %   q(3) => Angular speed of rotary arm q(1)_dot
        %   q(4) => Angular speed of pendulum arm q(2)_dot
    
        for i=1:n.trajs_training
            % Initial state
            q(1,1,i) = pi/2*(2*rand() - 1); % [-pi/2;pi/2]
            region_id = randi([2, n.regions-1]);   % [1;16]
            q(2,1,i) = -pi + angle_region *(region_id + rand() - 1);
            q(3,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
            q(4,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
            
            j = 1;
            while j <= n.steps
                q(:,j+1,i) = f_discrete(q(:,j,i), u(:,j,i));
                % if (abs(q(1,j+1,i)) >= pi || abs(q(2,j+1,i) - q(2,1,i)) >= 2*pi || abs(q(2,j+1,i)) >= pi)
                %     u(:,:,i) = 2*umax*rand(n.inputs, n.steps, 1) - umax;
                %     j = 1;
                %     retries = retries + 1;
                % else
                %     j = j+1;
                % end
                j = j+1;
            end
        end
        data.training.q = q;
        data.training.u = u;
        disp('Training data collected');

        % 2. TESTING DATA
        n.trajs_testing = n.trajs - n.trajs_training;
        angle_region = 2*pi/n.regions;
        q = zeros(n.states, n.steps+1, n.trajs_testing);
        u = 2*umax*rand(n.inputs, n.steps, n.trajs_testing) - umax; % [-umax, umax]
     
        for i=1:n.trajs_testing
            % Initial state
            q(1,1,i) = pi/2*(2*rand() - 1); % [-pi/2;pi/2]
            region_id = randi([2, n.regions-1]);   % [1;16]
            q(2,1,i) = -pi + angle_region *(region_id + rand() - 1);
            q(3,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
            q(4,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
            
            j = 1;
            while j <= n.steps
                q(:,j+1,i) = f_discrete(q(:,j,i), u(:,j,i));
                % if (abs(q(1,j+1,i)) >= pi || abs(q(2,j+1,i) - q(2,1,i)) >= 2*pi || abs(q(2,j+1,i)) >= pi)
                %     u(:,:,i) = 2*umax*rand(n.inputs, n.steps, 1) - umax;
                %     j = 1;
                %     retries = retries + 1;
                % else
                %     j = j+1;
                % end
                j = j+1;
            end
        end
        data.testing.q = q;
        data.testing.u = u;
        disp('Testing data collected');

    elseif strcmp(mode,'FRF')
        u = umax*prbs(8,8);
        q = zeros(n.states, size(u,1), n.trajs);
    
        % Reminder : States q :
        %   q(1) => Angle of rotary arm
        %   q(2) => Angle of pendulum arm
        %   q(3) => Angular speed of rotary arm q(1)_dot
        %   q(4) => Angular speed of pendulum arm q(2)_dot
    
        for i=1:n.trajs
            % Initial state
            q(1,1,i) = 0; %pi/2*(2*rand() - 1); % [-pi/2;pi/2]
            % region_id = randi(n.regions);   % [1;16]
            q(2,1,i) = 0; %-pi + (region_id - 0.5) * angle_region + (rand() - 0.5) * angle_region;
            q(3,1,i) = 0; %5*pi*(2*rand() - 1); % [-5*pi;5*pi]
            q(4,1,i) = 0; %5*pi*(2*rand() - 1); % [-5*pi;5*pi]
            
            j = 1;
            while j <= size(u,1)
                q(:,j+1,i) = f_discrete(q(:,j,i), u(j,1));
                j = j+1;
            end
        end
        data.q = q;
        data.u = u;
    else
        disp('This mode does not exist.');
    end
end
