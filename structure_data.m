function [data, n] = structure_data(tdms_data_path, n)
    data_tdms = tdmsread(tdms_data_path);
    n.trajs = numel(data_tdms);
    steps_plus1 = height(data_tdms{1});
    steps = steps_plus1 - 1;

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
        u = u';
        u = u(:,1:end-1);

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

    data.training.q = zeros(n.states, steps_plus1, n.trajs_training);
    data.training.u = zeros(n.inputs, steps,        n.trajs_training);
    data.testing.q  = zeros(n.states, steps_plus1,  n.trajs_testing);
    data.testing.u  = zeros(n.inputs, steps,         n.trajs_testing);

    for i = 1:n.trajs
        if i <= n.trajs_training
            data.training.q(:,:,i) = valid_trajs{i}.q;
            data.training.u(:,:,i) = valid_trajs{i}.u;
        else
            data.testing.q(:,:,i-n.trajs_training) = valid_trajs{i}.q;
            data.testing.u(:,:,i-n.trajs_testing)  = valid_trajs{i}.u;
        end
    end

    n.steps = steps;
    fprintf('Valid trajectories kept: %d / %d\n', n.trajs, numel(data_tdms));
end