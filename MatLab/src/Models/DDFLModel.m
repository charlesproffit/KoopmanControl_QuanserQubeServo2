classdef DDFLModel < handle

    properties (Access=public)

        nz
        ncontrolled

        Ts
        
        % Linear system
        A_c
        B_c

        f_lifting

        T
        gamma
        etta
        u

        v_h
        v_etta
        v_gamma
   
    end

    methods (Access=public)

        function obj = DDFLModel(f_lifting,nz,ncontrolled,Ts)
        
            obj.f_lifting = @(q)f_lifting(q);
        
            obj.nz = nz;
            obj.ncontrolled = ncontrolled;

            obj.Ts = Ts;
        
        end

        function DDFLIdentification(obj,model)

            phi = @(q)(obj.f_lifting(q));
            A = model.A_CT;
            B = model.B_CT;

            obj.v_h = zeros(obj.nz, 1);
            obj.v_h(2) = 1;

            M = [obj.v_h'; obj.v_h' * A];
            obj.T = @(q) M*phi(q);
            obj.v_gamma = obj.v_h'*(A^(obj.ncontrolled-1))*B;
            obj.gamma = @(q) obj.v_gamma*phi(q);
            obj.v_etta = -obj.v_h'*(A^(obj.ncontrolled));
            obj.etta = @(q) obj.v_etta*phi(q);

            obj.u = @(v,q)(v+obj.etta(q))/obj.gamma(q);

            obj.A_c = [
                [zeros(obj.ncontrolled-1,1), eye(obj.ncontrolled-1)];
                [0, zeros(1,obj.ncontrolled-1)];
            ];
        
            obj.B_c = [
                zeros(obj.ncontrolled-1,1);
                1;
            ];

        end

        % Continuous dynamics
        function xdot = continuousDynamics(obj,x,u)
            xdot = obj.A_c*x + obj.B_c*u;
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

    end

end