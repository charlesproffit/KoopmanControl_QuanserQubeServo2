function [FRF, data_lifted] = compute_frf(data, f, Ts, n, mode)
    if strcmp(mode,'FRF')
        n_used = n.lifted_states;
    else
        n_used = n.states;
    end
    % Step 1 : Compute lifted states
    lifted_q = zeros(n_used, size(data.u,1), n.trajs);
    for i = 1:n.trajs
        for j = 1:size(data.u,1)
            lifted_q(:,j,i) = f(data.q(:,j,i));
        end
    end
    data_lifted.q = lifted_q;
    data_lifted.u = data.u;

    t = 0:0.01:(size(data.u,1)-1)*0.01;
    figure;
    hold on
    plot(t, data_lifted.q(2,:,1));
    % plot(t, data_lifted.u);
    title("Pendulum arm angle");
    % legend("Response", "Input");
    hold off
    
    % Step 2 : Compute FRF between u and observables
    ufft = zeros(8,255,n.lifted_states);
    yfft = zeros(8,255,n.lifted_states);
    Gperiodic = zeros(n.lifted_states,1,255);
    % Gperiodic = zeros(1,255,n.lifted_states);
    ws = 0:(2*pi/(Ts*255)):(2*pi/(Ts*255))*254;
    
    for i = 1:n.lifted_states
        for j = 4:8
            ufft(j,:,i) = fft(data_lifted.u((j-1)*255+1:255*j));
            y = data_lifted.q(i,:,1);
            yfft(j,:,i) = fft(y((j-1)*255+1:255*j));
        end
    end

    for i = 1:n.lifted_states
        for j = 4:8
            Gperiodic(i,:,:) = Gperiodic(i,:,:) + yfft(j,:,i)./ufft(j,:,i);
        end
        Gperiodic(i,:,:) = Gperiodic(i,:,:)/5;
    end

    sysG = frd(Gperiodic,ws);

    FRF = sysG;
end