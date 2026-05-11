function compare_FL_models(M_MBFL, M_DDFL)
    %% 1. Vary q2, fix others
    q2_vec = linspace(-pi, pi, 300);
    q_base = [0; 0; 0; 0];

    gamma_ddfl  = zeros(1, length(q2_vec));
    etta_ddfl = zeros(1, length(q2_vec));
    gamma_mbfl    = zeros(1, length(q2_vec));
    etta_mbfl   = zeros(1, length(q2_vec));

    for k = 1:length(q2_vec)
        q = q_base;
        q(2) = q2_vec(k);
        gamma_ddfl(k)  = M_DDFL.gamma(q);
        etta_ddfl(k) = M_DDFL.etta(q);
        gamma_mbfl(k) = M_MBFL.gamma(q);
        etta_mbfl(k)   = M_MBFL.etta(q);
    end

    figure('Name','PFL model comparison: sweep over q2');

    subplot(1,3,1);
    plot(q2_vec*180/pi, gamma_mbfl,   'b-', 'LineWidth', 2); hold on;
    plot(q2_vec*180/pi, gamma_ddfl, 'r--','LineWidth', 2);
    xline( 90,'k:'); xline(-90,'k:');
    xlabel('q2 [deg]'); ylabel('Value');
    title('gamma');
    legend('Model-based','DDFL'); grid on;

    subplot(1,3,2);
    plot(q2_vec*180/pi, etta_mbfl,   'b-', 'LineWidth', 2); hold on;
    plot(q2_vec*180/pi, etta_ddfl, 'r--','LineWidth', 2);
    xline( 90,'k:'); xline(-90,'k:');
    xlabel('q2 [deg]'); ylabel('Value');
    title('etta');
    legend('Model-based','DDFL'); grid on;

    subplot(1,3,3);
    plot(q2_vec*180/pi, abs(gamma_mbfl - gamma_ddfl),   'b-', 'LineWidth',1.5); hold on;
    plot(q2_vec*180/pi, abs(etta_mbfl - etta_ddfl), 'r--','LineWidth',1.5);
    xline( 90,'k:'); xline(-90,'k:');
    xlabel('q2 [deg]'); ylabel('Absolute error');
    title('Absolute errors');
    legend('gamma error','etta error'); grid on;
    sgtitle('DDFL vs Model-based: sweep over q2 (q1=q3=q4=0)');

    %% 2. Sweep over q3 (q1dot) with q2=0 fixed
    q3_vec = linspace(-5*pi, 5*pi, 300);

    beta_ddfl2  = zeros(1, length(q3_vec));
    alpha_ddfl2 = zeros(1, length(q3_vec));
    gamma_mbfl2    = zeros(1, length(q3_vec));
    etta_mbfl2   = zeros(1, length(q3_vec));

    for k = 1:length(q3_vec)
        q = q_base;
        q(3) = q3_vec(k);
        beta_ddfl2(k)  = M_DDFL.gamma(q);
        alpha_ddfl2(k) = M_DDFL.etta(q);
        gamma_mbfl2(k) = M_MBFL.gamma(q);
        etta_mbfl2(k)  = M_MBFL.etta(q);
    end

    figure('Name','PFL model comparison: sweep over q3');

    subplot(1,3,1);
    plot(q3_vec, gamma_mbfl2,   'b-', 'LineWidth', 2); hold on;
    plot(q3_vec, beta_ddfl2, 'r--','LineWidth', 2);
    xlabel('q3=q1dot [rad/s]'); ylabel('Value');
    title('gamma');
    legend('Model-based','DDFL'); grid on;

    subplot(1,3,2);
    plot(q3_vec, etta_mbfl2,   'b-', 'LineWidth', 2); hold on;
    plot(q3_vec, alpha_ddfl2, 'r--','LineWidth', 2);
    xlabel('q3=q1dot [rad/s]'); ylabel('Value');
    title('etta');
    legend('Model-based','DDFL'); grid on;

    subplot(1,3,3);
    plot(q3_vec, abs(gamma_mbfl2 - beta_ddfl2),   'b-', 'LineWidth',1.5); hold on;
    plot(q3_vec, abs(etta_mbfl2 - alpha_ddfl2), 'r--','LineWidth',1.5);
    xlabel('q3=q1dot [rad/s]'); ylabel('Absolute error');
    title('Absolute errors');
    legend('gamma error','etta error'); grid on;
    sgtitle('DDFL vs Model-based: sweep over q3 (q1=q2=q4=0)');