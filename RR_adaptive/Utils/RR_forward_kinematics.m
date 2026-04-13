function [xi, xi_dot] = RR_forward_kinematics(q, dq, par_kin)
%RR_FORWARD_KINEMATICS end effector position and velocity

T_0_i = get_T_0_i(q, par_kin);
xi = T_0_i(1:2,end,end);

J = get_J(q, par_kin);
xi_dot = J(1:2,:) * dq;

end

