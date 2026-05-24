function M_DDFL_CT = compute_DDFL(f_lifting, M_BILINEAR_CT, n, tol)
    %% 1. System
    phi = @(q)(f_lifting(q));
    A_phiphi = M_BILINEAR_CT.A;
    B_phiphi = M_BILINEAR_CT.N;

    %% Test
    A = M_BILINEAR_CT.A;
    B = M_BILINEAR_CT.N;

    % For rho=2: C = [B, A*B]
    C_mat2 = [B];
    
    [UC, SC, VC] = svd(C_mat2');
    
    % Check each candidate v_h
    fprintf('%-5s %-15s %-15s\n', 'i', '||v_h^T * B||', '||v_h^T * A*B||');
    for i = size(UC,2):-1:1
        v_h_cand = UC(:,i);
        r1 = norm(v_h_cand' * B, 1);
        r2 = norm(v_h_cand' * A * B, 1);
        fprintf('%-5d %-15.6f %-15.6f\n', i, r1, r2);
    end

    %% 2. Find v
    v_h = zeros(n.lifted_states, 1);
    v_h(2) = 1;

    %% 3. Compute T, gamma and etta
    M = [v_h'; v_h' * A_phiphi];
    % M = zeros(2, n.lifted_states);
    % M(1,2) = 1; % To extract q1
    % M(2,4) = 1; % To extract q3
    T = @(q) M*phi(q);
    v_gamma = v_h'*(A_phiphi^(n.controlled_states-1))*B_phiphi;
    gamma = @(q) v_gamma*phi(q);
    v_etta = -v_h'*(A_phiphi^(n.controlled_states));
    etta = @(q) v_etta*phi(q);
    
    %% 4. Compute the dynamics DDFL model
    A_c = [
        [zeros(n.controlled_states-1,1), eye(n.controlled_states-1)];
        [0, zeros(1,n.controlled_states-1)];
    ];

    B_c = [
        zeros(n.controlled_states-1,1);
        1;
    ];

    M_DDFL_CT.A = A_c;
    M_DDFL_CT.B = B_c;
    M_DDFL_CT.T = @(q) T(q);
    M_DDFL_CT.gamma = @(q) gamma(q);
    M_DDFL_CT.etta = @(q) etta(q);
    M_DDFL_CT.u = @(q,v) (v + etta(q))/gamma(q);
    M_DDFL_CT.Tdot = @(q,v) (A_c*T(q) + B_c*v);

    Labview_weights.M_T = M;            % such that T(q) = M_T*phi(q)
    Labview_weights.v_gamma = v_gamma;  % such that gamma(q) = v_gamma*phi(q)
    Labview_weights.v_etta = v_etta;    % such that etta(q) = v_etta*phi(q)
    Labview_weights.v_h = v_h;
    save("data\Labview_weights.mat", "Labview_weights");
end