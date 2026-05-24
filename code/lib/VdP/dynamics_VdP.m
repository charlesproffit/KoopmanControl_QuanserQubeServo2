function [f_continuous, f_discrete] = dynamics_VdP(Ts, discretization_method)

    q1_dot =    @(q)(q(2));
    q2_dot =    @(q,u)(-q(1)+0.5*(1-q(1)^2)*q(2)+(1-q(2)^2)*u);
    
    % q_dot(t) = f_continuous(q(t), u(t))
    f_continuous = @(q,u)[q1_dot(q) ; q2_dot(q,u)];
    
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

end