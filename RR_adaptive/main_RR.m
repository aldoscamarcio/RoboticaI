clear all; close all; clc;
addpath('./functions');
addpath('./Utils');

%% RR Manipulator Definition
% RR params
g = 9.81;
a1 = 1.0;
a2 = 1.0;
a3 = 1.0;
m1 = 1.0;
m2 = 1.0;
m3 = 1.0;
l1 = a1/2;
l2 = a2/2;
l3 = a3/2;
Izz1 = (1/12)*m1*a1^2;
Izz2 = (1/12)*m2*a2^2;
Izz3 = (1/12)*m3*a3^2;
d1 = 0.5;
d2 = 0.5;
d3 = 0.5;
par_real = [m1;m2;m3; l1;l2;l3; Izz1;Izz2;Izz3; d1;d2;d3; a1;a2;a3;g];

% Model linear on parameters
% p1 = m;
% p2 = l*m;
% p3 = i + m*(l^2);
% p4 = d;
% % ----- INVERSE -----
% m = p1;
% l = p2/p1;
% i = p3 - (p2^2)/p1;
% d = p4;
p_real = [m1;m2;m3; l1*m1;l2*m2;l3*m3; Izz1+m1*(l1^2);Izz2+m2*(l2^2);Izz3+m3*(l3^2)];
par_kin = [a1, a2, a3];
par_nl = [a1, a2, a3, g];

% uncertainties
rng(1);
uncert = 0.15;
u_min = 1 - uncert;
u_max = 1 + uncert;
U = diag(u_min + rand(size(p_real))*(u_max - u_min));
% uncertain parameters
p_0 = U * p_real;

% Initialize jointstate
q0 = [pi/4; -pi/4; pi/4];
dq0 = [0.0 ; 0.0 ; 0.0];

%% Trajectory
% Lissajous trajectory parameters
time_gain = 0.5;                    % Riferimento velocemente variabile time_gain = 1
par_traj.xi0 = [0.8; 1.8];
par_traj.A = 0.5;
par_traj.B = 0.5;
a0 = 1;
b0 = 2; 
par_traj.a = time_gain * a0;
par_traj.b = time_gain * b0;
par_traj.d = 0;

% For Kinematic Inversion
q0_d = [pi/4; -pi/4; pi/4];
dq0_d = [0.0; 0.0; 0.0];
