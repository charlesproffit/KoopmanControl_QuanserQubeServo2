function trajs = compare_models_FL_VdP(M_MBFL, M_DDFL, f_discrete, n, Ts)
    umax = 0.5;
    noise = 2*umax*rand(n.inputs,n.steps,n.trajs_control) - umax;
    q_lti = zeros(n.states, n.steps, n.trajs_control);
    z_lti = zeros(n.states, n.steps, n.trajs_control);
    u_lti = zeros(n.inputs, n.steps, n.trajs_control);
    q_nl = zeros(n.states, n.steps, n.trajs_control);
    z_nl = zeros(n.states, n.steps, n.trajs_control);
    u_nl = zeros(n.inputs, n.steps, n.trajs_control);

    for i=1:n.trajs_control
        % Initial state
        q_nl(1,1,i) = 1*(2*rand() - 1);
        q_lti(1,1,i) = q_nl(1,1,i);
        q_nl(2,1,i) = 1*(2*rand() - 1);
        q_lti(2,1,i) = q_nl(2,1,i);
        z0 = M_DDFL.T(q_nl(:,1,i));
        z_lti(:,1,i) = z0;
        z_nl(:,1,i) = z0;

        j = 1;
        while j <= (n.steps - 1)
            u_lti(:,j,i) = noise(1,j,i);
            z_lti(:,j+1,i) = z_lti(:,j,i) + Ts*M_DDFL.Tdot(z_lti(:,j,i), u_lti(:,j,i));
        
            ett = M_DDFL.etta(q_nl(:,j,i));
            gamm = M_DDFL.gamma(q_nl(:,j,i));
            u_nl(:,j,i) = (u_lti(:,j,i) + ett)/gamm;
            q_nl(:,j+1,i) = f_discrete(q_nl(:,j,i), u_nl(:,j,i));
            z_nl(:,j+1,i) = M_DDFL.T(q_nl(:,j+1,i));
            j=j+1;
        end
    end
    trajs.q_nl = q_nl;
    trajs.z_nl = z_nl;
    trajs.u_nl = u_nl;
    trajs.q_lti = q_lti;
    trajs.z_lti = z_lti;
    trajs.u_lti = u_lti;
    disp('Trajectories with controller collected');
