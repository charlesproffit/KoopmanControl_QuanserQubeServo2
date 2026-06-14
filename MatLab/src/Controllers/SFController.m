classdef SFController < handle

    properties (Access=public)
        
        NLModel
        K
        
    end

    methods (Access=public)

        function obj = SFController(NLModel)

            obj.NLModel = NLModel;

        end

        function [u,v] = compute_u(obj, x,x_ref)

            u = -obj.K*(x-x_ref);
            v = 0;

        end

        % Design LQR
        function K = designLQR(obj,xeq,ueq,Q,R)

            [Ad,Bd] = obj.NLModel.linearizeDiscrete(xeq,ueq);

            K = dlqr(Ad,Bd,Q,R);

            obj.K = K;

        end

    end

end