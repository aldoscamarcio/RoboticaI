function handle=plot_robot(q,par_kin,real)

    T_0_i = get_T_0_i(q, par_kin);
    p0 = T_0_i(1:2,end,1);
    p1 = T_0_i(1:2,end,2);
    p2 = T_0_i(1:2,end,3);
    p3 = T_0_i(1:2, end, 4);  % End-effector (nuovo per 3DOF)
    % [p0, p1, p2] = RR_chain(q, params);

    if real
        plot(p0(1),p0(2),'ko','linewidth',4);
        % handle=plot([p0(1) p1(1) p2(1)],[p0(2) p1(2) p2(2)],'-ko','linewidth',2);
         handle = plot([p0(1) p1(1) p2(1) p3(1)], ...
                     [p0(2) p1(2) p2(2) p3(2)], ...
                     '-ro', 'LineWidth', 2, 'MarkerSize', 6);
    else
        plot(p0(1),p0(2),'bo','linewidth',4);
        % handle=plot([p0(1) p1(1) p2(1)],[p0(2) p1(2) p2(2)],'-bo','linewidth',2);
         handle = plot([p0(1) p1(1) p2(1) p3(1)], ...
                     [p0(2) p1(2) p2(2) p3(2)], ...
                     '--bo', 'LineWidth', 1.5, 'MarkerSize', 5);
    end
end