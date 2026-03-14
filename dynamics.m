function [f_continuous, f_discrete] = dynamics(Ts)

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
    l = L_p/2;
    b_p = 5*10^(-5);
    J_p = m_p*L_p^2/3;
    
    % States q :
    %   q(1) => Angle of rotary arm
    %   q(2) => Angle of pendulum arm
    %   q(3) => Angular speed of rotary arm q(1)_dot
    %   q(4) => Angular speed of pendulum arm q(2)_dot
    
    f1 =        @(q)(J_r+J_p*sin(q(2))^2);
    f2 =        @(q)(m_p*l*r*cos(q(2)));
    f3 =        @(q)(2*J_p*sin(q(2))*cos(q(2))*q(3)*q(4) - m_p*l*r*sin(q(2))*q(4)^2);
    f4 =        @(q)(-m_p*l*r*cos(q(2))/J_p);
    f5 =        @(q)(sin(q(2))*cos(q(2))*q(3)^2 - m_p*g*l*sin(q(2))/J_p - b_p*q(4)/J_p);
    f6 =        @(q,u)(k_m*(u-k_m*q(3))/R_m - b_r*q(3));
    
    q3_dot =    @(q,u)((f6(q,u) - f3(q) - f2(q)*f5(q))/(f1(q) + f2(q)*f4(q)));
    q4_dot =    @(q,u)(f4(q)*q3_dot(q,u) + f5(q));
    
    % q_dot(t) = f_continuous(q(t), u(t))
    f_continuous = @(q,u)[q(3) ; q(4) ; q3_dot(q,u) ; q4_dot(q,u)];
    
    % q_plus(t+Ts) = f_discrete(q(t), u(t))
    k1 =        @(q,u)[f_continuous(q,u)];
    k2 =        @(q,u)[f_continuous(q+Ts/2*k1(q,u),u)];
    k3 =        @(q,u)[f_continuous(q+Ts/2*k2(q,u),u)];
    k4 =        @(q,u)[f_continuous(q+Ts*k3(q,u),u)];
    f_discrete = @(q,u)(q + Ts/6*(k1(q,u) + 2*k2(q,u) + 2*k3(q,u) + k4(q,u)));

end