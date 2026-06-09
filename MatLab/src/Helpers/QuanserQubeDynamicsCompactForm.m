function [Mq, H, Phi, B] = QuanserQubeDynamicsCompactForm()
% States q :
%   q(1) => Angle of rotary arm
%   q(2) => Angle of pendulum arm
%   q(3) => Angular speed of rotary arm q(1)_dot
%   q(4) => Angular speed of pendulum arm q(2)_dot

    %% Parameters

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

    %% Compact matrix form

    m_11 = @(q)J_r+J_p*sin(q(2))^2;
    m_12 = @(q)m_p*l*r*cos(q(2)); %goes to 0 if pi is close to pi/2
    m_21 = @(q)m_12(q); %goes to 0 if pi is close to pi/2
    m_22 = @(q)J_p;
    Mq = @(q)[
        m_11(q) m_12(q);
        m_21(q) m_22(q)
    ]; 

    c_11 = @(q)2*J_p*sin(q(2))*cos(q(2))*q(4) + b_r + k_m^2/R_m;
    c_12 = @(q)-m_p*l*r*sin(q(2))*q(4);
    c_21 = @(q)-J_p*sin(q(2))*cos(q(2))*q(3);
    c_22 = @(q)b_p;

    phi_1 = @(q)0;
    phi_2 = @(q)m_p*g*l*sin(q(2));
    Phi = @(q)[
        phi_1(q);
        phi_2(q)
    ];

    b1 = k_m/R_m;
    b2 = 0;
    B = @(q)[
        b1; 
        b2
    ];

    h1 = @ (q) (c_11(q)*q(3)+c_12(q)*q(4));
    h2 = @ (q) (c_21(q)*q(3)+c_22(q)*q(4));
    H = @(q)[
        h1(q);
        h2(q)
    ];

end