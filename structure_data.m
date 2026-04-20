function [data, n] = structure_data(data_path, n, steps_length)
    if endsWith(data_path, ".mat")
        load(data_path, 'data');
        data_tdms = data;
        clear data;
    else
        data_tdms = tdmsread(data_path);
    end
    n.trajs = numel(data_tdms);

    % First pass: load all valid trajectories
    valid_trajs = {};
    valid_count = 0;
    for i = 1:n.trajs
        T = data_tdms{i};
        q1 = T{:,"q1"};
        q2 = T{:,"q2"};
        q3 = T{:,"q3"};
        q4 = T{:,"q4"};
        u  = T{:,"u"};

        q = [q1 q2 q3 q4]';
        q = q(:, 3:end);        % remove first zero row due to reinitializing the arrays
        u = u';   
        u = u(:, 2:end-2);      % remove first zero row and last input we dont need it
        
        if ~strcmp(steps_length, 'all')
            % Truncate to desired length
            q = q(:, 1:steps_length+1);
            u = u(:, 1:steps_length);
        end

        n.steps = size(u,2);

        % we discard trajectory if q2 leaves [-pi, pi]
        if any(abs(q(2,:)) > pi)
            continue;
        end
        
        valid_count = valid_count + 1;
        valid_trajs{valid_count} = struct('q', q, 'u', u);
    end

    n.trajs = numel(valid_trajs);
    n.trajs_training = fix(n.train_test_ratio * n.trajs);
    n.trajs_testing  = n.trajs - n.trajs_training;

    data.training.q = zeros(n.states, n.steps+1, n.trajs_training);
    data.training.u = zeros(n.inputs, n.steps,        n.trajs_training);
    data.testing.q  = zeros(n.states, n.steps+1,  n.trajs_testing);
    data.testing.u  = zeros(n.inputs, n.steps,         n.trajs_testing);

    % ar = 1:n.trajs;
    ar = randperm(n.trajs);

    for i = 1:n.trajs
        if i <= n.trajs_training
            data.training.q(:,:,i) = valid_trajs{ar(i)}.q;
            data.training.u(:,:,i) = valid_trajs{ar(i)}.u;
        else
            data.testing.q(:,:,i-n.trajs_training) = valid_trajs{ar(i)}.q;
            data.testing.u(:,:,i-n.trajs_training)  = valid_trajs{ar(i)}.u;
        end
    end

    fprintf('Valid trajectories kept: %d / %d\n', n.trajs, numel(data_tdms));
end