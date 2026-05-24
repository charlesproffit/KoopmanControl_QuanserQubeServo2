function M_MBFL = compute_MBFL(K_LQR)
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

    M_MBFL.A = [0, 1;0, 0];
    M_MBFL.B = [0;1];
    M_MBFL.etta = @(q)Alpha(q)/Beta(q);
    M_MBFL.gamma = @(q) 1/Beta(q);
    M_MBFL.T = @(q) Transformation(q);