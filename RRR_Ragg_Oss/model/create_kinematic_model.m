%% ------------------------------------------------------------------------------------ %%
%   author: Giorgio Simonini
%   date:   21/11/2023
%   info: This file creates the kinematic robot model
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
syms a1 a2 a3 real
par_kin = [a1, a2, a3];
p = sym('p%d_%d', [N_joints, 2], 'real');
disp('done!')

%% FRAMES CHAIN
% - Frames - %
frames_tr = sym([0,0,0; a1,0,0; a2,0,0; a3,0,0]);
frames_or = sym([0,0,0; 0,0,0; 0,0,0; 0,0,0]);
% - Frame chain - %
disp('Building frames chain ...')
T_0_i = sym(zeros(4,4,N_frames));
for index = 1 : N_frames    % joints + end-effector
	if index <= N_joints
		% joints frames
        R_tmp = R_x(frames_or(index,1)) * R_y(frames_or(index,2)) * R_z(frames_or(index,3) + q(index));
		d_tmp = frames_tr(index, :)';
		if index == 1
			T_0_i(:,:,index) = [R_tmp, d_tmp; 0,0,0, 1];
		else
			T_0_i(:,:,index) = T_0_i(:,:,index-1)*[R_tmp, d_tmp; 0,0,0, 1];
		end
	else
		% end-effector frame
		R_tmp = R_x(frames_or(index,1)) * R_y(frames_or(index,2)) * R_z(frames_or(index,3)); % there is no q
		d_tmp = frames_tr(index, :)';
		T_0_i(:,:,index) = T_0_i(:,:,index-1)*[R_tmp, d_tmp; 0,0,0, 1];
	end
end
T_0_i = simplify(T_0_i);
disp('done!')
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(T_0_i,'file','..\functions\get_T_0_i', 'Vars', {q, par_kin});
	disp('done!')
end

%% JACOBIAN
disp('Computing Jacobian ...')
% The position jacobian is the differentiation of the position
J_or = sym(zeros(3,N_joints));
ee_pos = T_0_i(1:3, 4, end);
J_pos = jacobian(ee_pos, q);
for index = 1 : N_joints
	J_or(:,index) = T_0_i(1:3, 3, index);
end
J = simplify([J_pos; J_or]);
disp('done!')
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(J,'file','..\functions\get_J', 'Vars', {q, par_kin});
	disp('done!')
end

%% JACOBIAN DERIVATIVE
disp('Computing Jacobian derivative ...')
J_dot = J;
for index = 1:size(J,2)
    J_dot(:,index) = jacobian(J(:,index),q)*dq;
end
disp('done!')
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(J_dot,'file','..\functions\get_J_dot', 'Vars', {q, dq, par_kin});
	disp('done!')
end

%% KINEMATIC REGRESSOR
xi = J * dq;
Y_kin = jacobian(xi, par_kin);
if CREATE_FUNCTIONS 
    disp(' Creating function ...')
    matlabFunction(Y_kin,'file','..\functions\get_Y_kin', 'Vars', {q, dq, par_kin});
	disp('done!')
end

%% SAVE MODEL
if SAVE_MODEL
    disp(' Saving model ...')
    model_kin.N_joints = N_joints;
	model_kin.N_frames = N_joints + 1;
	model_kin.T_0_i = T_0_i;
    model_kin.J = J;
	model_kin.J_dot = J_dot;
    model_kin.Y_kin = Y_kin;
    save ('model_kin.mat', 'model_kin');
    disp('done!')
end

%% OTHER
function R = R_x(phi)
    R = [1, 0, 0; 0, cos(phi), -sin(phi); 0, sin(phi), cos(phi)];
end

function R = R_y(theta)
    R = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
end

function R = R_z(psi)
    R = [cos(psi), -sin(psi), 0; sin(psi), cos(psi), 0; 0, 0, 1];
end
