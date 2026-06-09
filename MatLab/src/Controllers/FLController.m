classdef FLController < handle

    properties (Access=public)
        
        FLModel
        K
        
    end

    methods (Access=public)

        function obj = FLController(FLModel)

            obj.FLModel = FLModel;

        end

        function u = compute_u(obj, x,x_ref)
            
            err = obj.FLModel.T(x) - obj.FLModel.T(x_ref);
            v = -obj.K*err;
            u = obj.FLModel.u(v,x);

        end

        % Design LQR
        function K = designLQR(obj,Q,R)

            K = lqr(obj.FLModel.A_c,obj.FLModel.B_c,Q,R);

            obj.K = K;

        end

    end


end