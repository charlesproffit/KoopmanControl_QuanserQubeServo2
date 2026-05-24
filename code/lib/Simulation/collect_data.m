function [data,n] = collect_data(f_discrete, n)
    % 1. TRAINING DATA
    angle_region = 2*pi/n.regions;
    q = zeros(n.states, n.steps, n.trajs_training);
    noise = 2*n.umax*rand(n.inputs,n.steps,n.trajs_training) - n.umax;

    retries = 0;

    % Reminder : States q :
    %   q(1) => Angle of rotary arm
    %   q(2) => Angle of pendulum arm
    %   q(3) => Angular speed of rotary arm q(1)_dot
    %   q(4) => Angular speed of pendulum arm q(2)_dot

    for i=1:n.trajs_training
        % Initial state
        q(1,1,i) = 0;
        region_id = randi([2, n.regions-1]);    % [1;16]
        q(2,1,i) = -pi + angle_region *(region_id + rand() - 1);
        q(3,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
        q(4,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
        
        if i==1 
            q(1,1,1) = 0;
            q(2,1,1) = 0.8*pi;
            q(3,1,1) = 0;
            q(4,1,1) = 0;
        end

        j = 1;
        while j <= (n.steps - 1)
            q(:,j+1,i) = f_discrete(q(:,j,i), noise(:,j,i));
            % if (abs(q(1,j+1,i)) >= pi/2 || abs(q(2,j+1,i) - q(2,1,i)) >= 2*pi || abs(q(2,j+1,i)) >= pi)
            %     % 1. Re randomize initial states
            %     q(1,1,i) = 0;
            %     region_id = randi([2, n.regions-1]);    % [1;16]
            %     q(2,1,i) = -pi + angle_region *(region_id + rand() - 1);
            %     q(3,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
            %     q(4,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
            %     % 2. Re randomize input (ie randomize noise again)
            %     noise(:,:,i) = 2*n.umax*rand(n.inputs, n.steps, 1) - n.umax;
            %     j = 1;
            %     retries = retries + 1;
            % else
            %     j = j+1;
            % end
            j=j+1;
        end
    end
    data.training.q = q;
    data.training.u = noise; % ATTENTION ON A MIS LE CONTROLEUR, MTN ON L'OUBLIE !! => DONC DE noise A q !!
    disp('Training data collected');

    % 2. TESTING DATA
    angle_region = 2*pi/n.regions;
    q = zeros(n.states, n.steps, n.trajs_testing);
    noise = 2*n.umax*rand(n.inputs,n.steps,n.trajs_testing) - n.umax;
 
    for i=1:n.trajs_testing
        % Initial state
        q(1,1,i) = 0;
        region_id = randi([2, n.regions-1]);    % [1;16]
        q(2,1,i) = -pi + angle_region *(region_id + rand() - 1);
        q(3,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
        q(4,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
        
        if i==1 
            q(1,1,1) = 0;
            q(2,1,1) = 0.8*pi;
            q(3,1,1) = 0;
            q(4,1,1) = 0;
        end
         
        j = 1;
        while j <= (n.steps - 1)
            q(:,j+1,i) = f_discrete(q(:,j,i), noise(:,j,i));
            % if (abs(q(1,j+1,i)) >= pi/2 || abs(q(2,j+1,i) - q(2,1,i)) >= 2*pi || abs(q(2,j+1,i)) >= pi)
            %     % 1. Re randomize initial states
            %     q(1,1,i) = 0;
            %     region_id = randi([2, n.regions-1]);    % [1;16]
            %     q(2,1,i) = -pi + angle_region *(region_id + rand() - 1);
            %     q(3,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
            %     q(4,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
            %     % 2. Re randomize input
            %     noise(:,:,i) = 2*n.umax*rand(n.inputs, n.steps, 1) - n.umax;
            %     j = 1;
            %     retries = retries + 1;
            % else
            %     j = j+1;
            % end
            j=j+1;
        end
    end
    data.testing.q = q;
    data.testing.u = noise;
    disp('Testing data collected');
end
