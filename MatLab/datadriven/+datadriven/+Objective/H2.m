function [constraints, objective] = H2(plant, controller, omegas, opts)
%H2 Sets up a $H_2$ objective for a LFR plant
%   
%   This function sets up the $H_2$ objective for the generalised plant 
%   represented in Linear Fraction Representation (LFR) form. The functions
%   returns the LMI contraints and the objective as specified in [1].
% 
%   [constraints, objective] = H2(plant, controller, omegas)
% 
%   [constraints, objective] = H2(..., "epsNormalisation", eps) sets the eps for
%   normalisation of the constraint. Negative eps disables the normalisation.
% 
%   [constraints, objective] = H2(..., "epsSVD", eps) sets the eps for 
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

% Ensure that frequencies are ordered!
omegas = sort(omegas);

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

% ----- Integrator Setup -----
% integ = [W(1); diff(W(:))];
integ = [omegas(1); 0.5 * diff(omegas(:))] + [0.5 * diff(omegas(:)); 0]; % Trapezoidal

% ----- Frequency responses -----
Pf = freqresp(plant, omegas);
[Xf, Yf] = freqresp(controller, omegas);

% ----- SDPVAR definitions -----
norm_inf = norm(lft(frd(Pf, omegas, plant.Ts), ss(controller)), inf);

if controller.factorisation == "Right"
    Gamma = sdpvar(nz, nz, length(omegas), 'hermitian', 'complex');
    for i = 1:length(omegas)
        assign(Gamma(:,:,i), eye(nz, nz) * norm_inf^2);
    end
elseif controller.factorisation == "Left"
    Gamma = sdpvar(nw, nw, length(omegas), 'hermitian', 'complex');
    for i = 1:length(omegas)
        assign(Gamma(:,:,i), eye(nw, nw) * norm_inf^2);
    end
end

% ===== Constraints =====    
objective = 0;
constraints_tmp = cell(length(omegas), 1);

for i_W = 1:length(omegas)
    % Update objective
    objective = objective + integ(i_W) * trace(Gamma(:,:,i_W)) / pi;

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
    elseif controller.factorisation == "Left"
        % Normalisation of plant
        if (opts.epsSVD >= 0)
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
    end

    M = [Gamma(:, :, i_W) - Lambda, F;  F', PP];
    
    % Normalisation for numerical stability
    if (opts.epsNormalisation >= 0)
        [s,v,d] = svd(double(M));
        v = diag(v) + opts.epsNormalisation;
        L = s * diag(1./sqrt(v)) * d';
        M = ((L' * M * L) + (L' * M * L)') / 2;
    end

    constraints_tmp{i_W} = M >= 0;
end
constraints = [constraints_tmp{:}];

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