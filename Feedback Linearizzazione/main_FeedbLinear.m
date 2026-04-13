clear; clear all; clc;
syms theta1 theta2 theta3 q1 q2 q3 dq1 dq2 dq3 k1 k2 k3 dtheta1 dtheta2 dtheta3...
    v1 v2 v3 u1 u2 u3 real
syms eps1x eps1y eps2x eps2y eps3x eps3y x_des y_des x_dot_des y_dot_des...
     x_ddot_des y_ddot_des ddq1 ddq2 ddq3 b1 b2 b3 real

B = RR_mass_matrix(q1, q2, q3);
C = RR_coriolis_matrix(q1, q2, q3, dq1, dq2, dq3);
G = RR_gravity_vector(q1, q2, q3);

q = [q1; q2; q3];
dq = [dq1; dq2; dq3];
tau = [k1*(theta1 - q1); k2*(theta2 - q2); k3*(theta3 - q3)];
ing = [dtheta1; dtheta2; dtheta3];

q_dotdot = B \ (tau - C * dq - G);
q_dotdot = subs(q_dotdot, dtheta1, 0);
q_dotdot = subs(q_dotdot, dtheta2, 0);
q_dotdot = subs(q_dotdot, dtheta3, 0);

f = [dq1; dq2; dq3; q_dotdot(1); q_dotdot(2); q_dotdot(3); 0; 0; 0];

g = [zeros(3,3);
     zeros(3,3);
     eye(3)];

g1=g(:,1);
g2=g(:,2);
g3=g(:,3);

x = [q1; q2; q3; dq1; dq2; dq3; theta1; theta2; theta3]; % Vettore di Stato
v = [v1; v2; v3];  %Ingressi stabilizzanti

% Dinamica completa:
x_dot = f + g*ing;
%% Derivata uscita
a1 = 1; a2 = 1; a3 = 1;
xe = a1*cos(q1) + a2*cos(q1+q2) + a3*cos(q1+q2+q3);
ye = a1*sin(q1) + a2*sin(q1+q2) + a3*sin(q1+q2+q3);
y = [xe; ye; q3];


Lfh = jacobian(y, x)*f; % dy/dq Derivata prima
% Lg1h=jacobian(y,x)*g1;
% Lg2h=jacobian(y,x)*g2;
% Lg3h=jacobian(y,x)*g3;

L2fh = jacobian(Lfh,x)*f;  % Derivata seconda
% Lg1Lfh=jacobian(Lfh,x)*g1;
% Lg2Lfh=jacobian(Lfh,x)*g2;
% Lg3Lfh=jacobian(Lfh,x)*g3;

L3fh = jacobian(L2fh,x)*f; % Derivata terza
Lg1L2fh = jacobian(L2fh,x)*g1;
Lg2L2fh = jacobian(L2fh,x)*g2;
Lg3L2fh = jacobian(L2fh,x)*g3;

gam = (L3fh);

E = ([Lg1L2fh,Lg2L2fh,Lg3L2fh]);

matlabFunction(E, "File", "matE", "Vars", ...
                [q1, q2, q3, k1, k2, k3]);

matlabFunction(gam, 'File', 'matgam', ...
               'Vars', {q1, q2, q3, dq1, dq2, dq3, ...
                        theta1, theta2, theta3, ...
                        k1,k2,k3});


u = E \ ( - gam + v);

% Gain Ka Kv Kp
p = [-15; -10; -8];
A = [0 1 0; 0 0 1; 0 0 0];
b = [0; 0; 1];
K = acker(A, b, p);


%% Tentativo di verifica simbolica

Y = simplify( gam + E * u);
 
% Y =
% 
% v1
% v2
% v3