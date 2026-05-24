function M_DDFL_CT = compute_DDFL_VdP(f_lifting, M_BILINEAR_CT, n, tol)
    
    %% 1. System
    phi = @(q)(f_lifting(q));
    A_phiphi = M_BILINEAR_CT.A;
    B_phiphi = M_BILINEAR_CT.N;

    %% 2. Find v, M_FBL, V_Eta, V_Gamma
    V_0_Null = null(B_phiphi',"rational"); 
    V_1_Null = B_phiphi'*A_phiphi'*null(B_phiphi',"rational");
    for i = 1:size(V_1_Null,2)
        if norm(V_1_Null(:,i),2)^2 >= 1e-6
            break;
        end
    end
    v_h = V_0_Null(:,i);
    M_FBL = [v_h'; v_h'*A_phiphi];
    
    V_Eta = -(v_h'*(A_phiphi^2))';
    V_Gamma = (v_h'*A_phiphi*B_phiphi)';

    gamma = @(q) V_Gamma' * phi(q);
    etta  = @(q) V_Eta'  * phi(q);
    T     = @(q) M_FBL    * phi(q);


    %% 3. Compute the dynamics DDFL model
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

    Labview_weights.M_T = M_FBL;                % such that T(q) = M_T*phi(q)
    Labview_weights.v_gamma = V_Gamma;          % such that gamma(q) = v_gamma*phi(q)
    Labview_weights.v_etta = V_Eta;             % such that etta(q) = v_etta*phi(q)
    save("Labview_weights.mat", "Labview_weights");
end