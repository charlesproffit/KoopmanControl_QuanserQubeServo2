function [data,nsteps,ntrajs_training,ntrajs_testing] = structure_data(data_path, nstates, ninputs, train_test_ratio, steps_length, discarding)

    raw = tdmsread(data_path);

    valid_trajs  = {};
    valid_count  = 0;
    has_ref_glob = false;
    has_v_glob   = false;

    for i = 1:numel(raw)
        T = raw{i};
        q = [T{:,"q1"}, T{:,"q2"}, T{:,"q3"}, T{:,"q4"}]';
        u = T{:,"u"}';

        % Remove first two samples (array reinitialization problem in LabView)
        q = q(:, 3:end);
        u = u(:, 3:end);

        % Optional fields
        has_ref = ismember('ref', T.Properties.VariableNames);
        has_v   = ismember('v',   T.Properties.VariableNames);
        has_ref_glob = has_ref_glob || has_ref;
        has_v_glob   = has_v_glob   || has_v;

        if has_ref
            ref = T{:,"ref"}';
            ref = ref(:, 3:end);
        end
        if has_v
            v = T{:,"v"}';
            v = v(:, 3:end);
        end

        % Truncate if needed
        if ~strcmp(steps_length, 'all')
            q = q(:, 1:steps_length);
            u = u(:, 1:steps_length);
            if has_ref; ref = ref(:, 1:steps_length); end
            if has_v;   v   = v(:,   1:steps_length); end
        end

        % Discard if pendulum leaves [-pi, pi]
        if discarding && any(abs(q(2,:)) > pi)
            continue;
        end

        valid_count = valid_count + 1;
        traj = struct('X', q, 'U', u);
        if has_ref; traj.Ref = ref; end
        if has_v;   traj.V   = v;   end
        valid_trajs{valid_count} = traj;
    end

    clear data;

    ntrajs          = numel(valid_trajs);
    nsteps          = size(valid_trajs{1}.X, 2);
    ntrajs_training = fix(train_test_ratio * ntrajs);
    ntrajs_testing  = ntrajs - ntrajs_training;

    data.training.X = zeros(nstates, nsteps, ntrajs_training);
    data.training.U = zeros(ninputs, nsteps, ntrajs_training);
    data.testing.X  = zeros(nstates, nsteps, ntrajs_testing);
    data.testing.U  = zeros(ninputs, nsteps, ntrajs_testing);
    data.training.Ref = zeros(nstates, nsteps, ntrajs_training);
    data.testing.Ref  = zeros(nstates, nsteps, ntrajs_testing);
    data.training.V = zeros(1, nsteps, ntrajs_training);
    data.testing.V  = zeros(1, nsteps, ntrajs_testing);

    ar = randperm(ntrajs);
    for i = 1:ntrajs
        traj = valid_trajs{ar(i)};
        if i <= ntrajs_training
            split = 'training'; idx = i;
        else
            split = 'testing';  idx = i - ntrajs_training;
        end
        data.(split).X(:,:,idx) = traj.X;
        data.(split).U(:,:,idx) = traj.U;
        if has_ref_glob; data.(split).Ref(1,:,idx) = traj.Ref; end
        if has_v_glob;   data.(split).V(:,:,idx)   = traj.V;   end
    end

    fprintf('Valid trajectories kept: %d / %d\n', ntrajs, numel(raw));
end