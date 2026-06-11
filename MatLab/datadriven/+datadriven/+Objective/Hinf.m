function [constraints, objective] = Hinf(plant, controller, omegas, opts)
%HINF Sets up a $H_\infty$ objective for a LFR plant
%   
%   This function sets up the $H_\infty$ objective for the generalised plant
%   represented in Linear Fraction Representation (LFR) form. The function
%   returns the LMI contraints and the objective as specified in [1].
% 
%   [constraints, objective] = Hinf(plant, controller, omegas)
% 
%   [constraints, objective] = Hinf(..., "epsNormalisation", eps) sets the eps 
%   for normalisation of the constraint. Negative eps disables the
%   normalisation.
% 
%   [constraints, objective] = Hinf(..., "epsSVD", eps) sets the eps for
%   normalisation of the plant. Negative eps disables the normalisation.
% 
%   ----------------------------------------------------------------------------
%   Outputs
%     * constraints : List of LMI constraints in YALMIP framework.
%     * objective   : Objective in YALMIP framework.
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
    opts.epsNormalisation (1, 1) double = 1e-4;
    opts.epsSVD (1, 1) double = 1e-6;    
end

nw = length(plant.InputGroup.U1);
nu = length(plant.InputGroup.U2);
nz = length(plant.OutputGroup.Y1);
ny = length(plant.OutputGroup.Y2);

% Check if base assumption is satisfied
if controller.factorisation == "Right"
    if ny > nw
        error("Base assumption is not satisfied!!!");
    end
elseif controller.factorisation == "Left"
    if nu > nz
        error("Base assumption is not satisfied!!!");
    end
end

% ----- Frequency responses -----
Pf = freqresp(plant, omegas);
[Xf, Yf] = freqresp(controller, omegas);

% ----- SDPVAR definitions -----
gamma = sdpvar(1);
assign(gamma, norm(lft(frd(Pf, omegas, plant.Ts), ss(controller)), 'inf')^2);
    
% ===== Constraints =====
objective = gamma;
constraints_tmp = cell(length(omegas), 1);

for i_W = 1:length(omegas)
    % Prepare different values
    % [X, Y] = controller.evalfr(omegas(i_W));
    X = Xf(:, :, i_W);
    Y = Yf(:, :, i_W);

    P11 = Pf(plant.OutputGroup.Y1, plant.InputGroup.U1, i_W);
    P12 = Pf(plant.OutputGroup.Y1, plant.InputGroup.U2, i_W);
    P21 = Pf(plant.OutputGroup.Y2, plant.InputGroup.U1, i_W);
    P22 = Pf(plant.OutputGroup.Y2, plant.InputGroup.U2, i_W);
        
    if controller.factorisation == "Right"
        % Normalisation of plant
        if (opts.epsSVD >= 0)
            s = svds(P21, 1, "smallest") + opts.epsSVD; % Smallest singular value
            P21 = P21 / s;
            P12 = P12 * s;
        end
    
        Phi = pinv(P21) * (Y - P22*X);
        Phi_c = double(Phi);

        F = P11*Phi + P12*X;
        PP = Phi'*Phi_c + Phi_c'*Phi - Phi_c'*Phi_c;

        Psi = eye(nw) - pinv(P21)*P21;
        Lambda = (P11 * Psi) * (P11 * Psi)';
        Gamma = gamma * eye(nz);
        
    elseif controller.factorisation == "Left"
        if (opts.epsSVD >= 0)
            % Normalisation of plant
            s = svds(P12, 1, "smallest") + opts.epsSVD; % Smallest singular value
            P21 = P21 * s;
            P12 = P12 / s;
        end

        Phi = (Y - X*P22) * pinv(P12);
        Phi_c = double(Phi);

        F = (Phi*P11 + X*P21)';
        PP = Phi*Phi_c' + Phi_c*Phi' - Phi_c*Phi_c';

        Psi = eye(nz) - P12*pinv(P12);
        Lambda = (Psi * P11)' * (Psi * P11);
        Gamma = gamma * eye(nw);
    end
    
    M = [Gamma - Lambda, F;  F', PP];
    
    if (opts.epsNormalisation >= 0)
        % Normalisation for numerical stability
        [s,v,d] = svd(double(M));
        v = diag(v) + opts.epsNormalisation;
        L = s * sqrt(diag(1./v)) * d';
        M = ((L' * M * L) + (L' * M * L)') / 2;
    end
    
    constraints_tmp{i_W} = M >= eye(size(M))*1e-6;
end

constraints = [constraints_tmp{:}, gamma >= 0];

% Polygonal Chain
if controller.factorisation == "Right"
    if controller.ny == 1
        P22 = Pf(plant.OutputGroup.Y2, plant.InputGroup.U2, :);
        tmp = reshape(squeeze(Yf) - squeeze(Xf) .* squeeze(P22), 1, 1, []);
        constraints = [
            constraints;
            datadriven.Constraint.polygonal_chain( tmp );
            ];
    end
elseif controller.factorisation == "Left"
    if controller.nu == 1
        P22 = Pf(plant.OutputGroup.Y2, plant.InputGroup.U2, :);
        tmp = reshape(squeeze(Yf) - squeeze(P22) .* squeeze(Xf), 1, 1, []);
        constraints = [
            constraints;
            datadriven.Constraint.polygonal_chain( tmp );
            ];
    end
end

end