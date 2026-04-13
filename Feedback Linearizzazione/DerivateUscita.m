
syms q1(t) q2(t) q3(t) dq1(t) dq2(t) dq3(t) ddq1(t) ddq2(t) ddq3(t) theta1(t) theta2(t) theta3(t) dtheta1(t) dtheta2(t) dtheta3(t) real 
syms tau(t) k1 k2 k3 real

a1 = 1; a2 = 1; a3 = 1;
xe(t) = a1*cos(q1(t)) + a2*cos(q1(t)+q2(t)) + a3*cos(q1(t)+q2(t)+q3(t));
ye(t) = a1*sin(q1(t)) + a2*sin(q1(t)+q2(t)) + a3*sin(q1(t)+q2(t)+q3(t));
y(t) = [xe(t); ye(t); q3(t)];
Ke = [k1 , k2, k3];

x = [q1(t) q2(t) q3(t) dq1(t) dq2(t) dq3(t)  theta1(t) theta2(t) theta3(t)];
q(t) = [q1(t), q2(t), q3(t)]';
dq(t) = [dq1(t), dq2(t), dq3(t)]';
theta(t) = [theta1(t) theta2(t) theta3(t)]';
dtheta(t) = [dtheta1(t) dtheta2(t) dtheta3(t)]';

B = RR_mass_matrix(q1(t), q2(t), q3(t));
C_ = RR_coriolis_matrix(q1(t), q2(t), q3(t), dq1(t), dq2(t), dq3(t));
G_ = RR_gravity_vector(q1(t), q2(t), q3(t));

tau(t) = [k1*(theta1(t) - q1(t)); k2*(theta2(t) - q2(t)); k3*(theta3(t) - q3(t))];

ddq = B \ ((Ke * tau(t)) - C_ - G_);

ddq1(t)=ddq(1);
ddq2(t)=ddq(2);
ddq3(t)=ddq(3);


dy = diff(y(t),t);
dy = subs(dy, diff(q1(t), t) , dq1(t));
dy = subs(dy, diff(q2(t), t) , dq2(t));
dy = subs(dy, diff(q3(t), t) , dq3(t));

ddy = diff(dy,t);

ddy = subs(ddy, diff(q1(t), t) , dq1(t));
ddy = subs(ddy, diff(q2(t), t) , dq2(t));
ddy = subs(ddy, diff(q3(t), t) , dq3(t));
ddy = subs(ddy, diff(dq1(t), t) , ddq1(t));
ddy = subs(ddy, diff(dq2(t), t) , ddq2(t));
ddy = subs(ddy, diff(dq3(t), t) , ddq3(t));

ddy= simplify(ddy);

 


