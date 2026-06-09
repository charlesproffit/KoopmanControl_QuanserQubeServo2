classdef NLModel < handle

    properties(Access=public)

        % Dimensions
        nx
        nu
        nq

        % Sampling time
        Ts

        % Manipulator model
        Mq
        H
        Phi
        B_u

        % Dynamics
        f_nr
        g_nr
        f
        g

        % LQR data
        K_LQR = []
        x_eq = []
        u_eq = []

    end

    methods (Access=public)
        % Constructor
        function obj = NLModel(nx,nu,Ts,Mq,H,Phi,B_u)

            obj.nx = nx;
            obj.nu = nu;
            obj.nq = nx/2;

            obj.Ts = Ts;

            obj.Mq = Mq;
            obj.H = H;
            obj.Phi = Phi;
            obj.B_u = B_u;

            obj.f_nr = @(x)(-obj.Mq(x)\(obj.H(x)+obj.Phi(x)));
            obj.f = @(x)[
                x(obj.nq+1:end);
                obj.f_nr(x)
            ];

            obj.g_nr = @(x)(obj.Mq(x)\obj.B_u(x));
            obj.g = @(x)[
                zeros(obj.nq,obj.nu);
                obj.g_nr(x)
            ];

        end

        % Continuous dynamics
        function xdot = continuousDynamics(obj,x,u)

            if ~isempty(obj.K_LQR)
                u = u - obj.K_LQR*(x-obj.x_eq);
            end

            xdot = obj.f(x) + obj.g(x)*u;

        end

        % Euler step
        function xplus = stepEuler(obj,x,u)

            xplus = x + obj.Ts*obj.continuousDynamics(x,u);

        end

        % RK4 step
        function xplus = stepRK4(obj,x,u)

            k1 = obj.continuousDynamics(x,u);

            k2 = obj.continuousDynamics( x + obj.Ts/2*k1 , u);

            k3 = obj.continuousDynamics( x + obj.Ts/2*k2 , u);

            k4 = obj.continuousDynamics( x + obj.Ts*k3 , u);

            xplus = x + obj.Ts/6*( k1 + 2*k2 + 2*k3 + k4 );

        end

        % Linearization
        function [A,B] = linearize(obj,xeq,ueq)

            syms q1 q2 q3 q4 u real
            q = [q1; q2; q3; q4];
            
            %st M*q_dd+H+Phi=B*u
            f_syms = [
                    q3;
                    q4;
                    obj.f_nr(q)
                ];
            g_syms = [
                    0;
                    0;
                    obj.g_nr(q)
                ];
        
            x_dot = f_syms+g_syms*u;
        
            A = jacobian(x_dot, q);
            B = jacobian(x_dot, u);
        
            A_lin = simplify(subs(A, [q1 q2 q3 q4 u], [xeq.' ueq]));
            B_lin = simplify(subs(B, [q1 q2 q3 q4 u], [xeq.' ueq]));
        
            A = double(A_lin);
            B = double(B_lin);

        end

        % Discrete linearization
        function [Ad,Bd] = linearizeDiscrete(obj,xeq,ueq)

            [A,B] = obj.linearize(xeq,ueq);

            sysc = ss(A,B,eye(obj.nx), zeros(obj.nx,obj.nu));

            sysd = c2d(sysc,obj.Ts,'zoh');

            Ad = sysd.A;
            Bd = sysd.B;

        end

        % Design LQR
        function K = designLQR(obj,xeq,ueq,Q,R)

            [Ad,Bd] = obj.linearizeDiscrete(xeq,ueq);

            K = dlqr(Ad,Bd,Q,R);

            obj.K_LQR = K;
            obj.x_eq = xeq;
            obj.u_eq = ueq;

        end

    end

end