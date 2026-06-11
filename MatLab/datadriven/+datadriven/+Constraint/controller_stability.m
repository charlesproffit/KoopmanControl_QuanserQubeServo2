function [constraints] = controller_stability(controller, omegas, opts)
%controller_stability Ensure controller stablity
%   
%   constraints = controller_stability(controller, omegas)
% 
%   constraints = controller_stability(..., "epsNormalisation", eps) sets
%   eps for normalisation of the constraint. Negative eps disables the
%   normalisation.
%
%   ------------------------------------------------------------------------
%   Outputs
%     * constraints : List of LMI constraints in YALMIP framework.
% 
%   ------------------------------------------------------------------------
%   References
%   [1] P. Schuchert, V. Gupta, and A. Karimi, "Data-driven fixed-structure 
%       frequency-based $\mathcal{H}_2$ and $\mathcal{H}_\infty$ controller 
%       design," Automatica, vol. 160, p. 111398, Feb. 2024, 
%       doi: 10.1016/j.automatica.2023.111398.
% 
%   ------------------------------------------------------------------------
%   Copyright 2025 Vaibhav Gupta, DDMAC, EPFL (MIT License)
%

arguments (Input)
    controller datadriven.Controller.Controller
    omegas (:, 1) double
    opts.epsNormalisation (1, 1) double = 1e-4;
end

% Ensure that frequencies are ordered!
omegas = sort(omegas);

% ----- Frequency responses -----
[~, Yf] = freqresp(controller, omegas);
    
% ===== Constraints =====
constraints_tmp = cell(length(omegas), 1);

for i_W = 1:length(omegas)
    % Prepare different values
    % [~, Y] = controller.evalfr(omegas(i_W));
    Y = Yf(:, :, i_W);
    Y_c = double(Y);

    if controller.factorisation == "Right"
        M = Y'*Y_c + Y_c'*Y - Y_c'*Y_c; 
    elseif controller.factorisation == "Left"
        M = Y*Y_c' + Y_c*Y' - Y_c*Y_c';
    end

    if (opts.epsNormalisation >= 0)
        % Normalisation for numerical stability
        [s,v,d] = svd(double(M));
        v = diag(v) + opts.epsNormalisation;
        L = s * sqrt(diag(1./v)) * d';
        M = ((L' * M * L) + (L' * M * L)') / 2;
    end
    
    constraints_tmp{i_W} = M >= eye(size(M))*1e-6;
end
constraints = [constraints_tmp{:}];

if controller.nu == 1 && controller.ny == 1
    constraints = [
        constraints;
        datadriven.Constraint.polygonal_chain(Yf);
        ];
end

end