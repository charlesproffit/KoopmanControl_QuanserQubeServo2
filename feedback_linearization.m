function [M_DDFL_CT, M_DDFL_DT] = feedback_linearization(f_lifting, M_BILINEAR_CT, n, tol, Ts)
    % 1. Augment the system to have the system phi+ = A*phi+B*phi*u
    phi = @(q)[
        1;
        f_lifting(q);
    ];

    A = [
        [0, zeros(1,size(M_BILINEAR_CT.A,1))];
        [zeros(size(M_BILINEAR_CT.A,1),1), M_BILINEAR_CT.A];
    ];

    B = [
        [0, zeros(1,size(M_BILINEAR_CT.N,1))];
        [M_BILINEAR_CT.B, M_BILINEAR_CT.N];
    ];


    % 2. Find v
    Co = ctrb(A,B); % this constructs [B, A*B, ..., A^(n.lifted_states-1)*B]
    Co = Co(:,1:(n.states-1));
    [U, ~, ~] = svd(Co);
    v = zeros(size(U,1),1); % column vector
    for i=n.lifted_states:-1:1
        u = U(:,i); % this or this : U(i,:)' ??
        if norm(u'*A^(n.states-1)*B,1) >= tol
            v = u;
            break
        end
    end

    % 3. Compute T, gamma and etta
    % M = v' ;
    % for i=1:n.states-1
    %     M = [M; v'*A^(i)];
    % end
    % T = @(q) M*phi(q);
    M = zeros(4, n.lifted_states+1);
    M(1,2) = 1; % To extract q1
    M(2,4) = 1; % To extract q3
    T = @(q) M*phi(q);
    gamma = @(q) v'*A^(n.states-1)*B*phi(q);
    etta = @(q) -v'*A^(n.states)*phi(q);
    
    % 4. Compute the dynamics DDFL model
    A_c = [
        [zeros(n.states-1,1), eye(n.states-1)];
        [0, zeros(1,n.states-1)];
    ];

    B_c = [
        zeros(n.states-1,1);
        1;
    ];

    M_DDFL_CT.A = A_c;
    M_DDFL_CT.B = B_c;
    M_DDFL_CT.T = @(q) T(q);
    M_DDFL_CT.gamma = @(q) gamma(q);
    M_DDFL_CT.etta = @(q) etta(q);
    % M_DDFL_CT.u = @(q,v) (v + etta(q))/gamma(q);

    % Discretization
    sys_ct = ss(M_DDFL_CT.A, M_DDFL_CT.B, eye(size(M_DDFL_CT.A,1)), 0);
    sys_dt = c2d(sys_ct, Ts, 'zoh');
    M_DDFL_DT.A = sys_dt.A;
    M_DDFL_DT.B = sys_dt.B;
    M_DDFL_DT.T = @(q) T(q);
    M_DDFL_DT.gamma = @(q) gamma(q);
    M_DDFL_DT.etta = @(q) etta(q);
    % M_DDFL_DT.u = @(q,v) (v + etta(q))/gamma(q);

end