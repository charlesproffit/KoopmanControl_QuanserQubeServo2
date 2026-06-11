function [constraints, slack] = H2(plant, controller, omegas, opts)
%H2 Sets up a $H_2$ constraint for a LFR plant
%   
%   This function sets up the $H_2$ constraint for the generalised plant
%   represented in Linear Fraction Representation (LFR) form. The function 
%   returns the LMI contraints and the slack as specified in [1].
% 
%   [constraints, slack] = H2(plant, controller, omegas)
% 
%   [constraints, objective] = H2(..., "strict", false) solves the slacked 
%   version of the constraint.
% 
%   [constraints, slack] = H2(..., "epsNormalisation", eps) sets the eps for 
%   normalisation of the constraint. Negative eps disables the normalisation.
% 
%   [constraints, slack] = H2(..., "epsSVD", eps) sets the eps for normalisation
%   of the plant. Negative eps disables the normalisation.
% 
%   ----------------------------------------------------------------------------
%   Outputs
%     * constraints : List of LMI constraints in YALMIP framework.
%     * slack       : Slack in YALMIP framework or 0 if strict constraint
% 
%   ----------------------------------------------------------------------------
%   References
%   [1] P. Schuchert, V. Gupta, and A. Karimi, "Data-driven fixed-structure 
%       frequency-based $\mathcal{H}_2$ and $\mathcal{H}_\infty$ controller 
%       design," Automatica, vol. 160, p. 111398, Feb. 2024, 
%       doi: 10.1016/j.automatica.2023.111398.
% 
%   ----------------------------------------------------------------------------
%   Copyright 2025 Vaibhav Gupta, DDMAC, EPFL (MIT License)
%

arguments (Input)
    plant lti
    controller datadriven.Controller.Controller
    omegas (:, 1) double
    opts.strict (1, 1) logical = true
    opts.epsNormalisation (1, 1) double = 1e-4;
    opts.epsSVD (1, 1) double = 1e-6;  
end

[constraints, objective] = datadriven.Objective.H2(plant, controller, omegas, ...
    "epsNormalisation", opts.epsNormalisation, ...
    "epsSVD", opts.epsSVD);

if opts.strict
    slack = 0;
    constraints = [
        constraints;
        objective <= 1;
        ];
else
    slack = sdpvar(1);
    constraints = [
        constraints;
        objective <= 1 + slack;
        slack >= 0;
        ];
end

end