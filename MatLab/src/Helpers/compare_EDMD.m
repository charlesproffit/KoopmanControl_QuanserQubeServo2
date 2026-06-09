function comparison = compare_EDMD(dataset, X_BILIN, X_EDMD, X_LIN, X_BILIN_onestep, X_EDMD_onestep, X_LIN_onestep, BILIN, EDMD, LIN)

    Xreal = dataset.testing.X;
    U = dataset.testing.U;

    ntrajs = size(Xreal, 3);

    % Project lifted to original states
    X_BILIN_orig = reshape(BILIN.C * reshape(X_BILIN, BILIN.nz, []), BILIN.nx, [], ntrajs);
    X_EDMD_orig  = reshape(EDMD.C  * reshape(X_EDMD,  EDMD.nz,  []), EDMD.nx,  [], ntrajs);

    comparison.X_nonlinear = Xreal;
    comparison.X_BILIN     = X_BILIN_orig;
    comparison.X_EDMD      = X_EDMD_orig;
    comparison.X_LIN       = X_LIN;
    comparison.U           = U;
    comparison.U_training  = dataset.training.U;
    comparison.X_training  = dataset.training.X;

    % Testing errors in original state space
    err_BILIN = zeros(1, ntrajs);
    err_EDMD  = zeros(1, ntrajs);
    err_LIN   = zeros(1, ntrajs);
    for i = 1:ntrajs
        err_BILIN(i) = norm(X_BILIN_orig(:,:,i) - Xreal(:,:,i), 'fro');
        err_EDMD(i)  = norm(X_EDMD_orig(:,:,i)  - Xreal(:,:,i), 'fro');
        err_LIN(i)   = norm(X_LIN(:,:,i)        - Xreal(:,:,i), 'fro');
    end
    comparison.errors.BILIN_fro  = mean(err_BILIN);
    comparison.errors.EDMD_fro   = mean(err_EDMD);
    comparison.errors.LIN_fro    = mean(err_LIN);

    % Testing errors in lifted state space
    nsteps = size(Xreal, 2);
    err_BILIN_lifted = zeros(1, ntrajs);
    err_EDMD_lifted  = zeros(1, ntrajs);
    for i = 1:ntrajs
        lifted_true = zeros(BILIN.nz, nsteps);
        for t = 1:nsteps
            lifted_true(:,t) = BILIN.f_lifting(Xreal(:,t,i));
        end
        err_BILIN_lifted(i) = norm(X_BILIN(:,:,i) - lifted_true, 'fro');
        err_EDMD_lifted(i)  = norm(X_EDMD(:,:,i)  - lifted_true, 'fro');
    end
    comparison.errors.BILIN_fro_lifted = mean(err_BILIN_lifted);
    comparison.errors.EDMD_fro_lifted  = mean(err_EDMD_lifted);

    % Training errors from models
    comparison.errors.BILIN_fro_training        = BILIN.training_error_fro;
    comparison.errors.BILIN_fro_training_lifted = BILIN.training_error_fro_lifted;
    comparison.errors.EDMD_fro_training         = EDMD.training_error_fro;
    comparison.errors.EDMD_fro_training_lifted  = EDMD.training_error_fro_lifted;
    comparison.errors.LIN_fro_training          = LIN.training_error_fro;
    comparison.errors.LIN_fro_training_lifted   = LIN.training_error_fro_lifted;

    % SVDs and condition numbers
    comparison.svds.BILIN = BILIN.svds;
    comparison.svds.EDMD  = EDMD.svds;
    comparison.svds.LIN   = LIN.svds;
    comparison.cond.BILIN = BILIN.condition_number;
    comparison.cond.EDMD  = EDMD.condition_number;
    comparison.cond.LIN   = LIN.condition_number;

    
    % Project lifted to original states
    X_BILIN_orig_onstep = reshape(BILIN.C * reshape(X_BILIN_onestep, BILIN.nz, []), BILIN.nx, [], ntrajs);
    X_EDMD_orig_onstep  = reshape(EDMD.C  * reshape(X_EDMD_onestep,  EDMD.nz,  []), EDMD.nx,  [], ntrajs);

    comparison.X_BILIN_onestep     = X_BILIN_orig_onstep;
    comparison.X_EDMD_onestep      = X_EDMD_orig_onstep;
    comparison.X_LIN_onestep       = X_LIN_onestep;

    % Testing errors in original state space
    err_BILIN = zeros(1, ntrajs);
    err_EDMD  = zeros(1, ntrajs);
    err_LIN   = zeros(1, ntrajs);
    for i = 1:ntrajs
        err_BILIN(i) = norm(X_BILIN_orig_onstep(:,:,i) - Xreal(:,:,i), 'fro');
        err_EDMD(i)  = norm(X_EDMD_orig_onstep(:,:,i)  - Xreal(:,:,i), 'fro');
        err_LIN(i)   = norm(X_LIN_onestep(:,:,i)        - Xreal(:,:,i), 'fro');
    end
    comparison.errors.BILIN_fro_onestep = mean(err_BILIN);
    comparison.errors.EDMD_fro_onestep  = mean(err_EDMD);
    comparison.errors.LIN_fro_onestep   = mean(err_LIN);

    % Testing errors in lifted state space
    nsteps = size(Xreal, 2);
    err_BILIN_lifted = zeros(1, ntrajs);
    err_EDMD_lifted  = zeros(1, ntrajs);
    for i = 1:ntrajs
        lifted_true = zeros(BILIN.nz, nsteps);
        for t = 1:nsteps
            lifted_true(:,t) = BILIN.f_lifting(Xreal(:,t,i));
        end
        err_BILIN_lifted(i) = norm(X_BILIN_onestep(:,:,i) - lifted_true, 'fro');
        err_EDMD_lifted(i)  = norm(X_EDMD_onestep(:,:,i)  - lifted_true, 'fro');
    end
    comparison.errors.BILIN_fro_lifted_onestep = mean(err_BILIN_lifted);
    comparison.errors.EDMD_fro_lifted_onestep  = mean(err_EDMD_lifted);
end