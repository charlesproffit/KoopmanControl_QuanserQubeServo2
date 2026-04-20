function [M_DDFL] = feedback_linearization(f_lifting, M_BILINEAR, n, tol)
    % 1. Augment the system to have the system phi+ = A*phi+B*phi*u
    phi = @(q)[
        1;
        f_lifting(q);
    ];

    A = [
        [0, zeros(1,size(M_BILINEAR.A,1))];
        [zeros(size(M_BILINEAR.A,1),1), M_BILINEAR.A];
    ];

    B = [
        [0, zeros(1,size(M_BILINEAR.N,1))];
        [M_BILINEAR.B, M_BILINEAR.N];
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
    M = v' ;
    for i=1:n.states-1
        M = [M; v'*A^(i)];
    end
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

    M_DDFL.A = A_c;
    M_DDFL.B = B_c;
    M_DDFL.T = @(q) T(q);
    M_DDFL.gamma = @(q) gamma(q);
    M_DDFL.etta = @(q) etta(q);
    % M_DDFL.u = @(q,v) (v + etta(q))/gamma(q);

end