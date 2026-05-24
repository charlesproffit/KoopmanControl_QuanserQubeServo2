function [data,n] = collect_data_VdP(f_discrete, n, mode)
    if strcmp(mode, 'EDMD')
        % 1. TRAINING DATA
        q = zeros(n.states, n.steps, n.trajs_training);
        noise = 2*n.umax*rand(n.inputs,n.steps,n.trajs_training) - n.umax;%n.umax*randn(n.steps,n.inputs)'; % This format so same random values
    
        for i=1:n.trajs_training
            % Initial state
            q(1,1,i) = 0.9*(2*rand() - 1);%0.85;
            q(2,1,i) = 0.1*(2*rand() - 1);%-0.02;
            
            j = 1;
            while j <= (n.steps - 1)
                q(:,j+1,i) = f_discrete(q(:,j,i), noise(:,j,i));
                j=j+1;
            end
        end
        data.training.q = q;
        data.training.u = noise;
        disp('Training data collected');

        % 2. TESTING DATA
        q = zeros(n.states, n.steps, n.trajs_testing);
        noise = n.umax*randn(n.inputs,n.steps,n.trajs_testing);
     
        for i=1:n.trajs_testing
            % Initial state
            q(1,1,i) = 10*(2*rand() - 1);
            q(2,1,i) = 10*(2*rand() - 1);
             
            j = 1;
            while j <= (n.steps - 1)
                q(:,j+1,i) = f_discrete(q(:,j,i), noise(:,j,i));
                j=j+1;
            end
        end
        data.testing.q = q;
        data.testing.u = noise;
        disp('Testing data collected');
    end
end
