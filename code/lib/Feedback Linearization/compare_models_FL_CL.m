function compare_models_FL_CL(trajs_CL_MBFL, trajs_CL_DDFL, gen_title)
    
    figure;
    subplot(2,2,1); hold on;
    plot(trajs_CL_MBFL.q(1,:,1),'LineWidth',4,'Color','#ff0000');
    plot(trajs_CL_DDFL.q(1,:,1),':','LineWidth',4,'Color','#00ff00');
    title('q1');
    grid on; grid minor; box on;
    legend('Model FBL','Proposed Method','Location','southwest');
    
    subplot(2,2,2); hold on;
    plot(trajs_CL_MBFL.q(2,:,1),'LineWidth',4,'Color','#ff0000');
    plot(trajs_CL_DDFL.q(2,:,1),':','LineWidth',4,'Color','#00ff00');
    title('q2');
    grid on; grid minor; box on;
    
    subplot(2,2,3); hold on;
    plot(trajs_CL_MBFL.q(3,:,1),'LineWidth',4,'Color','#ff0000');
    plot(trajs_CL_DDFL.q(3,:,1),':','LineWidth',4,'Color','#00ff00');
    title('q3');
    grid on; grid minor; box on;
    
    subplot(2,2,4); hold on;
    plot(trajs_CL_MBFL.q(4,:,1),'LineWidth',4,'Color','#ff0000');
    plot(trajs_CL_DDFL.q(4,:,1),':','LineWidth',4,'Color','#00ff00');
    title('q4');
    grid on; grid minor; box on;
    sgtitle(gen_title);
end