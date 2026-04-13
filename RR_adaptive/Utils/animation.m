%% PLOTS

close all;
clear L;
f = figure;
f.WindowState = 'maximized';
pause(0.1)

% Getting Vectors from Sim
e_out = out.e_out.Data;
v_out = out.v_out.Data;
q_out = out.q_out.Data;
q_des_out = out.q_des_out.Data;
xi_des_out = out.xi_des_out.Data;
t_out = out.tout;

% post process
res = 0.1;
[t_sim, q_sim] = adjust_time(t_out,q_out,res);
[~, q_sim_des] = adjust_time(t_out,q_des_out,res);
[~, xi_des] = adjust_time(t_out,xi_des_out,res);

% errors
h(1) = subplot(2,2,1);
hold on
plot(t_out,e_out(1,:), 'LineWidth', 1);
plot(t_out,e_out(2,:), 'LineWidth', 1);
grid on
legend('$e_1$','$e_2$','Interpreter','latex')
title('Errors')

% lyapunov
h(1) = subplot(2,2,3);
hold on
plot(t_out,v_out, 'LineWidth', 1);
grid on
title('Lyapunov (V)')

% robot motion
h(2) = subplot(2,2,[2,4]);
hold on
axis equal
axis ([-2.0 3.0 -2.0 3.0])
grid on
title('Animation')

% initialize animation
t_h = text(2.2,2.2,['(' num2str(0) ')']);
t_xi = plot(xi_des(1,:,:),xi_des(2,:,:),'r');
L(1)=plot_robot(q_sim(:,1),par_kin,1);
L(2)=plot_robot(q_sim_des(:,1),par_kin,0);
legend(L, {'execution','reference'})

pause

for i=1:1:size(t_sim,2)
    
    delete(L(1))
    delete(L(2))
    delete(t_h)
    
    t_h = text(2.2,2.2,['(' num2str(t_sim(i)) ')']);
    hold on
    L(1)=plot_robot(q_sim(:,i),par_kin,1);
    L(2)=plot_robot(q_sim_des(:,i),par_kin,0);
    legend(L, {'execution','reference'})
    
    drawnow
    pause(0.01)
    
end