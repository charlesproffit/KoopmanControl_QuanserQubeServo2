classdef MBFLModel < handle
    properties (Access=public)
        ncontrolled
        Ts

        A_c
        B_c

        T
        gamma
        etta
        u
    end

    methods (Access=public)

        function obj = MBFLModel(K_LQR, Ts)
            obj.Ts = Ts;
            obj.ncontrolled = 2;

            mr = 0.095; mp = 0.024; r = 0.085; Lp = 0.129;
            br = 1e-3;  bp = 5e-5;
            Jr = mr*r^2/3; Jp = mp*Lp^2/3;
            km = 0.042; kt = km; Rm = 8.4;
            g_const = 9.81;
            l = Lp/2;
            n_q = 2;
            b1 = km/Rm;

            % Inertia / Coriolis / Gravity
            m_11  = @(q) Jr + Jp*sin(q(2))^2;
            m_12  = @(q) mp*l*r*cos(q(2));
            m_22  = @(q) Jp;
            c_11  = @(q,qd) 2*Jp*sin(q(2))*cos(q(2))*qd(2) + br + km^2/Rm;
            c_12  = @(q,qd) -mp*l*r*sin(q(2))*qd(2);
            c_21  = @(q,qd) -Jp*sin(q(2))*cos(q(2))*qd(1);
            c_22  = @(q,qd) bp;
            phi_2 = @(q) mp*g_const*l*sin(q(2));

            h1 = @(q,qd) c_11(q,qd)*qd(1) + c_12(q,qd)*qd(2);
            h2 = @(q,qd) c_21(q,qd)*qd(1) + c_22(q,qd)*qd(2);

            % Collocated partial FL
            m_11_bar = @(q) m_11(q) - m_12(q)/m_22(q)*m_12(q);
            h_1_bar = @(q,qd) h1(q,qd) - m_12(q)/m_22(q)*h2(q,qd);
            phi_1_bar = @(q) -m_12(q)/m_22(q)*phi_2(q);

            Alpha = @(x) (h_1_bar(x(1:n_q), x(n_q+1:2*n_q)) + phi_1_bar(x(1:n_q))) / b1 + K_LQR*x;
            Beta = @(x) m_11_bar(x(1:n_q)) / b1;

            obj.T = @(x) [x(1); x(3)];
            obj.etta = @(x) Alpha(x) / Beta(x);
            obj.gamma = @(x) 1/Beta(x);
            obj.u = @(v,x) (v+obj.etta(x))/obj.gamma(x);

            obj.A_c = [0, 1; 0, 0];
            obj.B_c = [0; 1];
        end

    end
end