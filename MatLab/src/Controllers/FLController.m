classdef FLController < handle

    properties (Access=public)
        
        mode = "control"

        FLModel
        K

        rand_v_amplitude
        
    end

    methods (Access=public)

        function obj = FLController(FLModel)

            obj.FLModel = FLModel;

        end

        function [u,v] = compute_u(obj, x,x_ref)
            
            err = obj.FLModel.T(x) - obj.FLModel.T(x_ref);
            v = -obj.K*err;
            if obj.mode ~= "control"
                v = v + obj.rand_v_amplitude*(2*rand()-1);
            end
            u = obj.FLModel.u(v,x);

        end

        % Design LQR
        function K = designLQR(obj,Q,R)

            K = lqr(obj.FLModel.A_c,obj.FLModel.B_c,Q,R);

            obj.K = K;

        end

        function switch_mode(obj, mode, amplitude)

            obj.mode = mode;
            obj.rand_v_amplitude = amplitude;

        end

    end


end