function [data, n] = structure_data(tdms_data_path, n)
    data_tdms = tdmsread(tdms_data_path);
    
    n.trajs = numel(data_tdms);
    
    steps_plus1 = height(data_tdms{1});
    steps = steps_plus1 - 1;
    
    n.trajs_training = fix(n.train_test_ratio*n.trajs);
    n.trajs_testing = n.trajs - n.trajs_training;

    data.training.q = zeros(n.states, steps_plus1, n.trajs_training);
    data.training.u = zeros(n.inputs, steps, n.trajs_training);

    data.testing.q = zeros(n.states, steps_plus1, n.trajs_testing);
    data.testing.u = zeros(n.inputs, steps, n.trajs_testing);

    for i = 1:n.trajs

        T = data_tdms{i};        % table
        M = table2array(T);      % convert to double
    
        % states
        q = M(:,1:4)';           % 4 × (steps+1)
    
        % input
        u = M(:,5)';             % 1 × (steps+1)
        u = u(:,2:end);          % remove first input
    
        if i <= n.trajs_training
            data.training.q(:,:,i) = q;
            data.training.u(:,:,i) = u;
        else
            data.testing.q(:,:,i-n.trajs_training) = q;
            data.testing.u(:,:,i-n.trajs_training) = u;
        end
        n.steps = size(data.training.u, 2);
    end
end