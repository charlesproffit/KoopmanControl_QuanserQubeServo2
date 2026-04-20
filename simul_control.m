function [trajs] = simul_control(f_discrete, umax, n, K_LQR, M_DDFL, K, z_ref)
    % 1. TRAINING DATA
    n.trajs_control = 100;
    angle_region = 2*pi/n.regions;
    q = zeros(n.states, n.steps+1, n.trajs_control);
    noise = 2*umax*rand(n.inputs,n.steps,n.trajs_control) - umax;
    u = zeros(n.inputs, n.steps, n.trajs_control);

    % Reminder : States q :
    %   q(1) => Angle of rotary arm
    %   q(2) => Angle of pendulum arm
    %   q(3) => Angular speed of rotary arm q(1)_dot
    %   q(4) => Angular speed of pendulum arm q(2)_dot

    for i=1:n.trajs_control
        % Initial state
        q(1,1,i) = 0; %pi/2*(2*rand() - 1); % [-pi/2;pi/2]
        region_id = randi([2, n.regions-1]);    % [1;16]
        q(2,1,i) = -pi + angle_region *(region_id + rand() - 1);
        q(3,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
        q(4,1,i) = 5*pi*(2*rand() - 1); % [-5*pi;5*pi]
        
        j = 1;
        while j <= n.steps
            z = M_DDFL.T(q(:,j,i));
            err = z - z_ref;
            v = - K*err;
            ett = M_DDFL.etta(q(:,j,i));
            gamm = M_DDFL.gamma(q(:,j,i));
            u(:,j,i) = (v + ett)/gamm - K_LQR*q(:,j,i);
            q(:,j+1,i) = f_discrete(q(:,j,i), u(:,j,i));
            j=j+1;
        end
    end
    trajs.q = q;
    trajs.u = u; % ATTENTION ON A MIS LE CONTROLEUR, MTN ON L'OUBLIE !! => DONC DE noise A q !!
    disp('Trajectories with controller collected');


end