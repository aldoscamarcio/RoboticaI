%% ------------------------------------------------------------------------------------ %%
%   author: Giorgio Simonini
%   date:   21/11/2023
%   info: This file creates the dynamic robot model
%   convention:
%       - Links start from 0, first bring the base frame to the first
%           joint, last bring the last joint to the end effector
%       - Frames start from 0, there are as many frames as links
%       - Frames have the origins on the phisical joints (0 and ee
%           excluded)
%       - Joints start from 1, frames and joints have same origins and
%           indexes
%       - Link Li brings the Joint {Ji} to the Joint {Ji+1} with a 
%           rototraslational matrix, defined by a traslation and a 
%           rotation using the XYZ parametrization
%       - Center of mass of link Li is expressed in frame {Si}
%       - Inertia tensor of link Li is defined on the center of mass in the
%           frame {Si}
%       - the joint variables are included in the trasformation matrices
% ---------------------------------------------------------------------------------------%

%% ----- INIT ----- %%
clear all
addpath('..\functions\')
% - defines - %
CREATE_FUNCTIONS = 1;
SAVE_MODEL = 1;
% - hyperparameters - %
N_joints = 3;               % number of joints
N_frames = N_joints + 1;    % joints + ee

%% VARIABLES AND PARAMETERS
disp('Variables and parameters initialization ...')
% - variables - %
q = sym('q%d', [N_joints,1], 'real');
dq = sym('dq%d', [N_joints,1], 'real');
ddq = sym('ddq%d', [N_joints,1], 'real');
% - parameters - %
syms a1 a2 a3 g real
par_nl = [a1, a2, a3 g];
masses = sym('m%d', [N_joints,1], 'real');              % Link Masses
center_of_masses = sym('l%d', [N_joints,1], 'real');    % Center of Masses
inertia = sym('i%d', [N_joints,1], 'real');             % Inertia Tensors
damping = sym('d%d', [N_joints,1], 'real');             % joint dampings
par = [masses, center_of_masses, inertia, damping];     % parameters set
par = reshape(par,[N_joints*size(par,2),1]);
par = [par; par_nl'];
p = sym('p%d_%d', [N_joints, 3], 'real');
disp('done!')

%% load kinematic model
load model_kin.mat
T_0_i = model_kin.T_0_i;

%% INERTIA MATRIX
disp('Computing Inertia Matrix ...')
% - center of mass jacobians - %
J_cm_i = sym(zeros(6, N_joints, N_joints));
J_i_or = sym(zeros(3,N_joints));
for index = 1 : N_joints
    l_tmp = [center_of_masses(index,1); 0; 0];
    d_0_i = T_0_i(1:3,4,index) + T_0_i(1:3,1:3,index)*l_tmp;
	d_0_i = simplify(d_0_i);
    J_i_pos = simplify(jacobian(d_0_i, q));
	for index2 = 1 : index
		if index <= N_joints
        	J_i_or(:,index2) = T_0_i(1:3, 3, index2);
		end
    end
    J_cm_i(:,:,index) = [J_i_pos; J_i_or];
end
% - inertia matrix - %
M = sym(zeros(N_joints));
for index = 1 : N_joints
    I_tmp = inertia(index,1);
    I_cm_i = [0, 0, 0; 0, 0, 0; 0, 0, I_tmp];
    I_cm_i = simplify(T_0_i(1:3,1:3,index) * I_cm_i * (T_0_i(1:3,1:3,index)'));
    M_cm_i = [masses(index)*eye(3), zeros(3); zeros(3), I_cm_i];
    M = M + J_cm_i(:,:,index)' * M_cm_i * J_cm_i(:,:,index);
end
disp('simplifying ...')
M = simplify(M);
disp('done!')
if CREATE_FUNCTIONS 
    disp(' Creating Inertia function ...')
    matlabFunction(M,'file','..\functions\get_M', 'Vars', {q, par});
    disp('done!')
end

%% CORIOLIS MATRIX
disp('Computing Coriolis Matrix ...')
C = sym(zeros(N_joints));
for index = 1 : 1 : N_joints
    for index2 = 1 : 1 : N_joints
        gamma_sum = sym(0);
        for index3 = 1 : 1 : N_joints
            gamma = 1/2 * (diff(M(index,index2), q(index3)) + diff(M(index,index3), q(index2)) - diff(M(index2,index3), q(index)));
            gamma = simplify(gamma);
            gamma_sum = ( gamma_sum + gamma ) * dq(index3);
        end
        clear gamma;
        C(index,index2) = gamma_sum;
    end
end
toc
disp('simplifying ...')
C = simplify(C);
disp('done!')
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(C,'file','..\functions\get_C', 'Vars', {q, dq, par});
    disp('done!')
end

%% LINK DAMPING MATRIX
disp('Computing Link damping Matrix ...')
D = diag(damping);
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(D,'file','..\functions\get_D', 'Vars', {par});
end
disp('done!')

%% GRAVITY TERMS
disp('Computing Gravity terms ...')
G = sym(zeros(N_joints,1));
grav = [0; -g; 0];
for index = 1 : N_joints
    G = G + J_cm_i(:,:,index)'*[grav*masses(index);zeros(3,1)];
end
G = simplify(G);
disp('done!')
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(G,'file','..\functions\get_G', 'Vars', {q, par});
    disp('done!')
end

%% MODEL LINEAR IN PARAMETERS
disp('Creating linear parameters ...')
% p1 = m;
% p2 = l*m;
% p3 = i + m*(l^2);
% % ----- INVERSE -----
% m = p1;
% l = p2/p1;
% i = p3 - (p2^2)/p1;
masses_p = p(:, 1);                         % Link Masses
center_of_masses_p = p(:,2)./p(:,1); 	    % Center of Masses
inertia_p = p(:,3) - (p(:,2).^2)./p(:,1);	% Inertia Tensor
% damping_p = p(:,4);                       % joint Dampings
p = reshape(p,[N_joints*3,1]);
disp('done!')

%% LINEAR INERTIA MATRIX
disp('computing linear Inertia Matrix ...')
% - center of mass jacobians - %
J_cm_i_p = sym(zeros(6, N_joints, N_joints));
J_i_or = sym(zeros(3,N_joints));
for index = 1 : N_joints
    l_tmp = [center_of_masses_p(index,1); 0; 0];
    d_0_i = T_0_i(1:3,4,index) + T_0_i(1:3,1:3,index)*l_tmp;
	d_0_i = simplify(d_0_i);
    J_i_pos = simplify(jacobian(d_0_i, q));
	for index2 = 1 : index
		if index <= N_joints
        	J_i_or(:,index2) = T_0_i(1:3, 3, index2);
        end
    end
    J_i_or = simplify(J_i_or);
    J_cm_i_p(:,:,index) = [J_i_pos; J_i_or];
end
Bp = sym(zeros(N_joints));
for index = 1 : N_joints
    I_tmp = inertia_p(index,1);
    I_cm_i = [0, 0, 0; 0, 0, 0; 0, 0, I_tmp];
    I_cm_i = simplify(T_0_i(1:3,1:3,index) * I_cm_i * (T_0_i(1:3,1:3,index)'));
    M_cm_i = [masses_p(index)*eye(3), zeros(3); zeros(3), I_cm_i];
    Bp = Bp + J_cm_i_p(:,:,index)' * M_cm_i * J_cm_i_p(:,:,index);
end
disp('simplifying ...')
Bp = simplify(Bp);
disp('done!')
if CREATE_FUNCTIONS 
    disp(' Creating Inertia function ...')
    matlabFunction(Bp,'file','..\functions\get_Mp', 'Vars', {q, p, par_nl});
    disp('done!')
end

%% LINEAR CORIOLIS MATRIX
disp('Computing linear Coriolis Matrix ...')
Cp = sym(zeros(N_joints));
for index = 1 : 1 : N_joints
    for index2 = 1 : 1 : N_joints
        gamma_sum = sym(0);
        for index3 = 1 : 1 : N_joints
            gamma = 1/2 * (diff(Bp(index,index2), q(index3)) + diff(Bp(index,index3), q(index2)) - diff(Bp(index2,index3), q(index)));
            gamma = simplify(gamma);
            gamma_sum = ( gamma_sum + gamma ) * dq(index3);
        end
        clear gamma;
        Cp(index,index2) = gamma_sum;
    end
end
disp('simplifying ...')
Cp = simplify(Cp);
disp('done!')
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Cp,'file','..\functions\get_Cp', 'Vars', {q, dq, p, par_nl});
    disp('done!')
end

%% LINEAR GRAVITY TERMS
disp('Computing linear Gravity terms ...')
Gp = sym(zeros(N_joints,1));
grav = [0; -g; 0];
for index = 1 : N_joints
    Gp = Gp + J_cm_i_p(:,:,index)'*[grav*masses_p(index);zeros(3,1)];
end
Gp = simplify(Gp);
disp('done!')
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Gp,'file','..\functions\get_Gp', 'Vars', {q, p, par_nl});
    disp('done!')
end

%% MODELS AND REGRESSOR
disp('Computing regressor ...')
dq_r = sym('dq_r%d', [N_joints,1], 'real');
ddq_r = sym('ddq_r%d', [N_joints,1], 'real');
% - MODEL - %
f = simplify(Bp*ddq_r + Cp*dq_r + Gp);
% - Regressor - %
Y = jacobian(f, p);
disp('done!')
disp('simplifying ...')
Y = simplify(Y);
Y = collect(expand(Y));     % p/p=1
disp('done!')
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Y,'file','..\functions\get_Y', 'Vars', {q, dq, dq_r, ddq_r, par_nl});
    disp('done!')
end

if SAVE_MODEL
    disp(' Saving model ...')
    model.N_joints = N_joints;
	model.N_frames = N_frames;
	model.par = par;
    model.par_nl = par_nl;
    model.p = p;
    model.M = M;
    model.Mp = Bp;
    model.C = C;
    model.Cp = Cp;
    model.D = D;
    model.G = G;
    model.Gp = Gp;
	model.T_0_i = model_kin.T_0_i;
    model.J = model_kin.J;
    model.J_dot = model_kin.J_dot;
    model.Y_kin = model_kin.Y_kin;
    model.Y = Y;
    save ('model.mat', 'model');
    disp('done!')
end
