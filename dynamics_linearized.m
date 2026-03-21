function [LIN_CT, LIN_DT] = dynamics_linearized(q_eq, u_eq, Ts)
    
    % Symbolic variables to allow for symbolic derivation
    syms q1 q2 q3 q4 u real
    q = [q1; q2; q3; q4];

    % General parameters
    g = 9.81;
    k_t = 0.042;
    k_m = 0.042;
    R_m = 8.4;
    
    % Rotary arm parameters
    m_r = 0.095;
    r = 0.085;
    b_r = 10^(-3);
    J_r = m_r*r^2/3;
    
    % Pendulum arm parameters
    m_p = 0.024;
    L_p = 0.129;
    l = 0.0645;
    b_p = 5*10^(-5);
    J_p = m_p*L_p^2/3;
    
    % States q :
    %   q(1) => Angle of rotary arm
    %   q(2) => Angle of pendulum arm
    %   q(3) => Angular speed of rotary arm q(1)_dot
    %   q(4) => Angular speed of pendulum arm q(2)_dot

    f1 = (J_r + J_p*sin(q2)^2);
    f2 = (m_p*l*r*cos(q2));
    f3 = (2*J_p*sin(q2)*cos(q2)*q3*q4 - m_p*l*r*sin(q2)*q4^2);
    f4 = (-m_p*l*r*cos(q2)/J_p);
    f5 = (sin(q2)*cos(q2)*q3^2 - m_p*g*l*sin(q2)/J_p - b_p*q4/J_p);
    f6 = (k_m*(u - k_m*q3)/R_m - b_r*q3);
    
    q3_dot = (f6 - f3 - f2*f5) / (f1 + f2*f4);
    q4_dot = f4*q3_dot + f5;
    
    f = [q3 ; q4 ; q3_dot ; q4_dot];

    A = jacobian(f, q);
    B = jacobian(f, u);

    A_lin = simplify(subs(A, [q1 q2 q3 q4 u], [q_eq.' u_eq]));
    B_lin = simplify(subs(B, [q1 q2 q3 q4 u], [q_eq.' u_eq]));

    A_num = double(A_lin);
    B_num = double(B_lin);

    LIN_CT.A = A_num;
    LIN_CT.B = B_num;

    % Discretization with ZOH method
    sys_ct = ss(LIN_CT.A, LIN_CT.B, [], []);
    sys_dt = c2d(sys_ct, Ts, 'zoh');
    LIN_DT.A = sys_dt.A;
    LIN_DT.B = sys_dt.B;
end