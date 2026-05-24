function [trajs] = simul_control(f_discrete, n, M_DDFL, K, q_ref)
    
    % Reference in transformed space
    z_ref = M_DDFL.T(q_ref);
    
    n.regions = 10;
    % 1. TRAINING DATA
    q = zeros(n.states, n.steps, n.trajs_control);
    u = zeros(n.inputs, n.steps, n.trajs_control);

    % Reminder : States q :
    %   q(1) => Angle of rotary arm
    %   q(2) => Angle of pendulum arm
    %   q(3) => Angular speed of rotary arm q(1)_dot
    %   q(4) => Angular speed of pendulum arm q(2)_dot
    
    angle_region = 2*pi/n.regions;
    for i=1:n.trajs_control
        % Initial state
        q(1,1,i) = 1/3*pi*(2*rand() - 1);
        q(2,1,i) = 2/3*pi*(2*rand() - 1);
        q(3,1,i) = 0;
        q(4,1,i) = 0;
        
        if i==1 
            q(1,1,1) = 0;
            q(2,1,1) = 0.8*pi;
            q(3,1,1) = 0;
            q(4,1,1) = 0;
        end
        
        j = 1;
        while j <= (n.steps - 1)
            z = M_DDFL.T(q(:,j,i));
            err = z - z_ref;
            v = - K*err;
            ett = M_DDFL.etta(q(:,j,i));
            gamm = M_DDFL.gamma(q(:,j,i));
            u(:,j,i) = (v + ett)/gamm;
            q(:,j+1,i) = f_discrete(q(:,j,i), u(:,j,i));
            j=j+1;
        end
    end
    trajs.q = q;
    trajs.u = u; % ATTENTION ON A MIS LE CONTROLEUR, MTN ON L'OUBLIE !! => DONC DE noise A q !!
    disp('Trajectories with controller collected');


end