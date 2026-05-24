function [f_continuous, f_discrete, f_discrete_with_LQR, parameters] = dynamics(Ts, discretization_method, K_LQR)

    % General parameters
    g = 9.81;
    k_t = 0.042;
    k_m = k_t;
    R_m = 8.4;
    
    % Rotary arm parameters
    m_r = 0.095;
    r = 0.085;
    b_r = 10^(-3);
    J_r = m_r*r^2/3; % around the pivot point
    
    % Pendulum arm parameters
    m_p = 0.024;
    L_p = 0.129;
    l = L_p/2;
    b_p = 5*10^(-5);
    J_p = m_p*L_p^2/3; % around the pivot point
    
    % States q :
    %   q(1) => Angle of rotary arm
    %   q(2) => Angle of pendulum arm
    %   q(3) => Angular speed of rotary arm q(1)_dot
    %   q(4) => Angular speed of pendulum arm q(2)_dot

    % f1 =        @(q)(J_r+J_p*sin(q(2))^2);
    % f2 =        @(q)(m_p*l*r*cos(q(2)));
    % f3 =        @(q)(2*J_p*sin(q(2))*cos(q(2))*q(3)*q(4) - m_p*l*r*sin(q(2))*q(4)^2);
    % f4 =        @(q)(-m_p*l*r*cos(q(2))/J_p);
    % f5 =        @(q)(sin(q(2))*cos(q(2))*q(3)^2 - m_p*g*l*sin(q(2))/J_p - b_p*q(4)/J_p);
    % f6 =        @(q,u)(k_m*(u-k_m*q(3))/R_m - b_r*q(3));
    % D  =        @(q)(f1(q) + f2(q)*f4(q));
    % 
    % q3_dot =    @(q,u)((f6(q,u) - f3(q) - f2(q)*f5(q))/D(q));
    % q4_dot =    @(q,u)(f4(q)*q3_dot(q,u) + f5(q));
    % 
    % q_dot(t) = f_continuous(q(t), u(t))
    % f_continuous = @(q,u)[q(3) ; q(4) ; q3_dot(q,u) ; q4_dot(q,u)];
    % 
    % q_plus(t+Ts) = f_discrete(q(t), u(t))
    % if strcmp(discretization_method, 'rk4')
    %     k1 =        @(q,u)[f_continuous(q,u)];
    %     k2 =        @(q,u)[f_continuous(q+Ts/2*k1(q,u),u)];
    %     k3 =        @(q,u)[f_continuous(q+Ts/2*k2(q,u),u)];
    %     k4 =        @(q,u)[f_continuous(q+Ts*k3(q,u),u)];
    %     f_discrete = @(q,u)(q + Ts/6*(k1(q,u) + 2*k2(q,u) + 2*k3(q,u) + k4(q,u)));
    % elseif strcmp(discretization_method, 'euler')
    %     f_discrete = @(q,u)(q+Ts*f_continuous(q,u));
    % end
    % 
    % parameters.f1 = @(q)f1(q);
    % parameters.f2 = @(q)f2(q);
    % parameters.f3 = @(q)f3(q);
    % parameters.f4 = @(q)f4(q);
    % parameters.f5 = @(q)f5(q);
    % parameters.f6 = @(q,u)f6(q,u);
    % parameters.D  = @(q)D(q);


    %% Qube System
    m_11 = @(q)J_r+J_p*sin(q(2))^2;
    m_12 = @(q)m_p*l*r*cos(q(2)); %goes to 0 if pi is close to pi/2
    m_21 = @(q)m_12(q); %goes to 0 if pi is close to pi/2
    m_22 = @(q)J_p;
    c_11 = @(q)2*J_p*sin(q(2))*cos(q(2))*q(4) + b_r + k_m^2/R_m;
    c_12 = @(q)-m_p*l*r*sin(q(2))*q(4);
    c_21 = @(q)-J_p*sin(q(2))*cos(q(2))*q(3);
    c_22 = @(q)b_p;
    phi_1 = @(q)0;
    phi_2 = @(q)m_p*g*l*sin(q(2));
    % n_u = 1; %number of inputs
    % n_q = 2; %nb of joint coordinates
    b1 = k_m/R_m;
    b2 = 0;
    h1 = @ (q) (c_11(q)*q(3)+c_12(q)*q(4));
    h2 = @ (q) (c_21(q)*q(3)+c_22(q)*q(4));
    Mq = @(q)[
        m_11(q) m_12(q);
        m_21(q) m_22(q)
    ];
    % Cq = @(q)[
    %     c_11(q) c_12(q);
    %     c_21(q) c_22(q)
    % ];
    Phi = @(q)[
        phi_1(q);
        phi_2(q)
    ];
    B = [b1; b2];
    H = @(q)[
        h1(q);
        h2(q)
    ];
    %st M*q_dd+H+Phi=B*u
    f = @(q)[
            q(3);
            q(4);
            -Mq(q)^(-1)*(H(q)+Phi(q))
        ];
    g = @(q)[
            0;
            0;
            Mq(q)^(-1)*B
        ];

    parameters.f = @(q)f(q);
    parameters.g = @(q)g(q);
    parameters.g_not_redundant = @(q)(Mq(q)^(-1)*B);
    parameters.f_not_redundant = @(q)(-Mq(q)^(-1)*(H(q)+Phi(q)) - parameters.g_not_redundant(q)*K_LQR*q);

    % q_dot(t) = f_continuous(q(t), u(t))
    f_continuous = @(q,u)(f(q)+g(q)*u);
    
    % q_plus(t+Ts) = f_discrete(q(t), u(t))
    if strcmp(discretization_method, 'rk4')
        k1 =        @(q,u)[f_continuous(q,u)];
        k2 =        @(q,u)[f_continuous(q+Ts/2*k1(q,u),u)];
        k3 =        @(q,u)[f_continuous(q+Ts/2*k2(q,u),u)];
        k4 =        @(q,u)[f_continuous(q+Ts*k3(q,u),u)];
        f_discrete = @(q,u)(q + Ts/6*(k1(q,u) + 2*k2(q,u) + 2*k3(q,u) + k4(q,u)));
    elseif strcmp(discretization_method, 'euler')
        f_discrete = @(q,u)(q+Ts*f_continuous(q,u));
    end
    f_discrete_with_LQR = @(q,u) f_discrete(q,u - K_LQR*q);


end