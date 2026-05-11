function M_MBFL = model_based_feedback_linearization(K_LQR)
    %% Physical parameters
    mr = 0.095; %[kg] rotary arm mass
    r = 0.085;  % [m] rotary arm length (8.5 cm)
    br = 1e-3;  % [N*m*s/rad] rotary arm damping
    Jr = mr*r^2/3; % [kg*m^2] rotary arm inertia
    mp = 0.024;    % [kg] pendulum mass (24 g)
    Lp = 0.129;    % [m] pendulum length (12.9 cm)
    bp = 5e-5;    % [N*m*s/rad] pendulum damping
    Jp = mp*Lp^2/3; % [kg*m^2] pendulum inertia about pivot
    km = 0.042;    % [N*m/A] Motor back-emf constant
    kt = km;    % [N*m/A] motor torque constant
    Rm = 8.4;    % [Ohm] motor resistance
    g_const = 9.81; % [m/s^2] gravitational acceleration
    l = Lp/2;

    %% Qube System
    m_11 = @(q)Jr+Jp*sin(q(2))^2;
    m_12 = @(q)mp*l*r*cos(q(2)); %goes to 0 if pi is close to pi/2
    m_21 = @(q)m_12(q); %goes to 0 if pi is close to pi/2
    m_22 = @(q)Jp;
    c_11 = @(q,qdot)2*Jp*sin(q(2))*cos(q(2))*qdot(2) + br + km^2/Rm;
    c_12 = @(q,qdot)-mp*l*r*sin(q(2))*qdot(2);
    c_21 = @(q,qdot)-Jp*sin(q(2))*cos(q(2))*qdot(1);
    c_22 = @(q,qdot)bp;
    phi_1 = @(q)0;
    phi_2 = @(q)mp*g_const*l*sin(q(2));
    % n_u = 1; %number of inputs
    n_q = 2; %nb of joint coordinates
    b1 = km/Rm;
    % b2 = 0;
    % inputs = [b1;b2];
    h1 = @ (q,qdot) (c_11(q,qdot)*qdot(1)+c_12(q,qdot)*qdot(2));
    h2 = @ (q,qdot) (c_21(q,qdot)*qdot(1)+c_22(q,qdot)*qdot(2));

    % For partial feedback linearization- collocated
    m_11_bar = @(q) m_11(q)- (m_12(q)*(1/m_22(q))*m_21(q));
    h_1_bar = @(q,qdot) h1(q,qdot)- m_12(q)*(1/m_22(q))*h2(q,qdot);
    phi_1_bar = @(q) phi_1(q)- m_12(q)*(1/m_22(q))*phi_2(q);
    Alpha = @(x)(h_1_bar(x(1:n_q),x(n_q+1:2*n_q))+phi_1_bar(x(1:n_q)) )/b1 + K_LQR*x;
    Beta = @(x) m_11_bar(x(1:n_q))/b1;
    Transformation = @(x) [x(1); x(3)]; % q1 and q1dot
    
    % For partial feedback linearization- non-collocated
    % m_12_tilde = @(q) m_12(q)- m_11(q)*m_22(q)/m_21(q);
    % h_1_tilde = @(q,qdot) h1(q,qdot)- m_11(q)*h2(q,qdot)/m_21(q);
    % phi_1_tilde = @(q) phi_1(q)- m_11(q)*phi_2(q)/m_21(q);
    % Alpha_NC = @(x) (h_1_tilde(x(1:n_q), x(n_q+1:2*n_q)) + phi_1_tilde(x(1:n_q)) )/b1;
    % Beta_NC = @(x) (m_12_tilde(x(1:n_q)) )/b1;
    % Transformation_NC = @(x) [x(2); x(4)]; % it’s q2 and q2dot

    % %to pass more easily to qube_dynamics.m
    % SYS.m11 = m_11; SYS.m12 = m_12; SYS.m21 = m_21; SYS.m22 = m_22;
    % SYS.c11 = c_11; SYS.c12 = c_12; SYS.c21 = c_21; SYS.c22 = c_22;
    % SYS.phi1 = phi_1; SYS.phi2 = phi_2;
    % SYS.b1 = b1;
    % 
    % q = x(1:2);
    % qd = x(3:4);
    % Mq = [SYS.m11(q) SYS.m12(q);
    % SYS.m21(q)
    % SYS.m22(q)];
    % Cq = [
    %     SYS.c11(q,qd) SYS.c12(q,qd);
    %     SYS.c21(q,qd) SYS.c22(q,qd)
    % ];
    % Phi = [
    %     SYS.phi1(q);
    %     SYS.phi2(q)
    % ];
    % B = [SYS.b1; 0];
    % qdd = Mq \ (B*u- Cq*qd- Phi);
    % xdot = [qd; qdd];

    M_MBFL.etta = @(q)Alpha(q)/Beta(q);
    M_MBFL.gamma = @(q) 1/Beta(q);
    M_MBFL.T = @(q) Transformation(q);