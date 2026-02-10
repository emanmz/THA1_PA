%% THA 1 - PA 
clear all; close all; clc;

%% 1a. rotational matrix SO3 -> equivalent axis angle representation

function [angle, axis] = rot2axisangle(R)
% check if so3
if isSO3(R) == 0
    error("Not a valid SO3 matrix");
end

if norm(R - eye(3), 'fro') < 1e-6 % CASE A: if R = I -> theta = 0
    error("Axis is undefined! Singularity!");
elseif trace(R)+1 < 1e-6 % CASE B: if tr(R) = -1 -> theta = pi
    angle = pi;
    % check what diagonal element is largest to avoid division by 0
    [~, large] = max(diag(R));
    if large == 1
        % w = 1/sqrt(2(1+r11) [ 1+ r11, r21, r31]
        n = 1/(sqrt(2*(1+R(1,1))));
        axis = n.*[R(1,1)+1; R(2,1); R(3,1)];
    elseif large == 2
        % w = 1/sqrt(2(1+r22) [ r12, 1+r22, r32]
        n = 1/(sqrt(2*(1+R(2,2))));
        axis = n.*[R(1,2); 1+ R(2,2); R(3,2)];
    else
        % w = 1/sqrt(2(1+r33) [ r13, r23, 1+r33]
        n = 1/(sqrt(2*(1+R(3,3))));
        axis = n.*[R(1,3); R(2,3); 1+ R(3,3)];
    end
else % CASE C: General Case 
    % theta = acos ( 1/2 ( trR - 1) E [0, pi)
    angle = acos(0.5*(trace(R)-1));
    % w = 1/sintheta ( R - R')
    axis = [R(3,2)-R(2,3);
        R(1,3)-R(3,1);
        R(2,1)-R(1,2)] / (2 * sin(angle));
end
end

%% 1b. rotational matrix -> Quaternion

function [q0, q1, q2, q3] = rot2quat(R)
% check if so3
if isSO3(R) == 0
    error("Not a valid SO3 matrix");
end
% any function from -pi to pi no singularities
% q0 = 1/2 sqrt( r11 + r22 + r33 + 1)
% q1 = 1/2 sgn(r32 - r23) sqrt( r11 - r22 - r33 + 1)
% q2 = 1/2 sgn(r13 - r31) sqrt( r22 - r33 - r11 + 1)
% q3 = 1/2 sgn(r21 - r12) sqrt( r33 - r11 - r22 + 1)

q0 = 0.5 * sqrt(R(1,1) + R(2,2) + R(3,3) + 1);
q1 = 0.5 * sgn(R(3,2) - R(2,3)) * sqrt(R(1,1) - R(2,2) - R(3,3) + 1);
q2 = 0.5 * sgn(R(1,3) - R(3,1)) * sqrt(R(2,2) - R(3,3) - R(1,1) + 1);
q3 = 0.5 * sgn(R(2,1) - R(1,2)) * sqrt(R(3,3) - R(1,1) - R(2,2) + 1);
end
 

%% 1c. rotation matrix -> ZYZ

function [phi, theta, psi] = rot2zyz(R)
if isSO3(R) == 0
    error("Not a valid SO3 matrix");
end

theta = atan2(sqrt(R(1,3)^2 + R(2,3)^2), R(3,3)); % in the range of 0, pi

% singularity sin(theta) = 0 or pi
if abs(sin(theta)) < 1e-6
    error("Singularity! Axes are parralel");
else
    % general case
    phi = atan2(R(2,3), R(1,3));
    psi = atan2(R(3,2), - R(3,1));

end
end

%% 1c. rotation matrix -> roll pitch yaw

function [roll, pitch, yaw] = rot2xyz(R)
if isSO3(R) == 0
    error("Not a valid SO3 matrix");
end 

% pitch = atan2(-r31, sqrt(r32^2 + r33^2))
pitch = atan2(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2));

if abs(cos(pitch)) < 1e-6 % When pitch is 90 or -90 deg singularity!
    error("Singularity!");
else 
    % roll = atan2(r32, r33)
    roll = atan2(R(3,2), R(3,3));
    % yaw = atan2(r21, r11)
    yaw = atan2(R(2,1), R(1,1));
end 
end 
%% 2a. axis angle -> rotation

function [R] = axisangle2rot(axis, angle)

% normalize axis so it's a unit vector
% || w || = 1

if norm(axis)> 0
    axis = axis / norm(axis);
end 

ux = axis(1);
uy = axis(2);
uz = axis(3);

% skew symmetric matrix  
K = [0, -uz, uy; uz, 0, -ux; -uy, ux, 0];

% rodriguez formula 
% R=I+(sintheta )K+(1costheta )K^2)
% w^2 = w w' - I (avoid matrix multiplication)
R = cos(angle) * eye(3) + (1 - cos(angle))* (axis * axis') + sin(angle) * K;

end 

%% 2b. quat -> rotation
function [R] = quat2rot(q0, q1, q2, q3)

% a unit quaternion
mag = sqrt(q0^2 + q1^2 + q2^2 + q3^2);
q0 = q0/mag; q1 = q1/mag; q2 = q2/mag; q3 = q3/mag;

% 2. Calculate matrix elements based on source formulas:
% Diagonal elements
r11 = q0^2 + q1^2 - q2^2 - q3^2; %
r22 = q0^2 - q1^2 + q2^2 - q3^2; %
r33 = q0^2 - q1^2 - q2^2 + q3^2; %

% Off-diagonal elements (Row 1)
r12 = 2*(q1*q2 - q0*q3); %
r13 = 2*(q0*q2 + q1*q3); %

% Off-diagonal elements (Row 2)
r21 = 2*(q0*q3 + q1*q2); %
r23 = 2*(q2*q3 - q0*q1); %

% Off-diagonal elements (Row 3)
r31 = 2*(q1*q3 - q0*q2); %
r32 = 2*(q0*q1 + q2*q3); %

% 3. Assemble the SO(3) Rotation Matrix
R = [r11, r12, r13;
    r21, r22, r23;
    r31, r32, r33];
end


%% Rotation matrix check

function SO3 = isSO3(R)

% 3 by 3 matrix
is3by3 = all(size(R) == [3 3]);
if ~is3by3
    SO3 = false;
    disp("NOT a Square 3 x 3 Rotational Matrix");
    return;
end 

% orthogonal: R'R = I
isOrthog = norm(R' * R - eye(3)) < 1e-6;

%  determinant: det(R) = 1
% If det(R) = -1, its an improper rotation (a reflection)
isDetOne = abs(det(R) - 1) < 1e-6;

SO3 = is3by3 && isOrthog && isDetOne;
end

%% sgn function 

function s = sgn(x)
if x >= 0 
    s = 1; 
else 
    s = -1;
end 
end 

%% Main / Testing functions for PA Q1 & Q2
clear all; close all; clc;

% Define test cases
test_names = {'90 deg X', '90 deg Y', '90 deg Z', '180 deg Z', 'Identity'};
matrices = {
    [1 0 0; 0 0 -1; 0 1 0], ... % 90 deg X
    [0 0 1; 0 1 0; -1 0 0], ... % 90 deg Y
    [0 -1 0; 1 0 0; 0 0 1], ... % 90 deg Z
    [-1 0 0; 0 -1 0; 0 0 1], ...% 180 deg Z (Case B Singularity)
    eye(3)                      % Identity (Case A Singularity)
};

for i = 1:length(matrices)
    fprintf('-----------Testing Case: %s -----------\n', test_names{i});
    R = matrices{i};
    
    % 1. Test Axis-Angle
    try
        disp("Rotational Matrix");disp(R);
        [ang, ax] = rot2axisangle(R);
        fprintf('Rotational 2 Axis Angle: Angle = %.4f, Axis = [%.2f, %.2f, %.2f]\n', ang, ax);
        disp("Builtin: Rotation Matrix 2 Axis Angle"); disp(rotm2axang(R));
        R_rec = axisangle2rot(ax, ang);
        disp("Axis Angle 2 Rotational"); disp(R_rec);
    catch ME
        fprintf('Axis-Angle Status: %s\n', ME.message);
    end
    
    % 2. Test Quaternions (Formula from)
    try
        [q0, q1, q2, q3] = rot2quat(R);
        fprintf('Quaternion: [%.4f, %.4f, %.4f, %.4f]\n', q0, q1, q2, q3);
        disp("Builtin: Rotation Matrix 2 Quat"); disp(rotm2quat(R));
        R_rec = quat2rot(q0, q1, q2, q3);
        disp("Quat 2 Rotational"); disp(R_rec);
    catch ME
        fprintf('Quaternion Status: %s\n', ME.message);
    end
    
    % 3. Test ZYZ (Formula from)
    try
        [phi, theta, psi] = rot2zyz(R);
        fprintf('ZYZ Euler: [%.4f, %.4f, %.4f]\n', phi, theta, psi);
        disp("Builtin: Rotation Matrix 2 ZYZ"); disp(rotm2eul(R, 'ZYZ'));
    catch ME
        fprintf('ZYZ Status: %s\n', ME.message);
    end
    
    % 4. Test XYZ (Roll-Pitch-Yaw)
    try
        [roll, pitch, yaw] = rot2xyz(R);
        fprintf('XYZ Euler: [%.4f, %.4f, %.4f]\n', roll, pitch, yaw);
        disp("Builtin: Rotation Matrix 2 XYZ"); disp(rotm2eul(R, 'XYZ'));
    catch ME
        fprintf('ZYX Status: %s\n', ME.message);
    end
    fprintf('\n');
end

%% Main for PA Q3 

% inputs 
T1 = [1 0 0 2; 0 1 0 0; 0 0 1 0; 0 0 0 1];
q = [0; 2; 0];
s_hat = [0; 0; 1];
h = 2; 

% screw axis matrix need toconvert point, direction and pitch to se3
S_mat = screw2Skew(q, s_hat, h);


% configurations 
thetas = [0, pi/4, pi/2, 3*pi/4, pi];
labels = {'0', 'pi/4', 'pi/2', '3pi/4', 'pi'};


% plot axis at all configs
figure; hold on; grid on; axis equal; view(3);
xlabel('X'); ylabel('Y'); zlabel('Z'); title('Rigid Body Motion Along Screw Axis');
T_configs = cell(1, 5);
for i = 1:5
    T_step = MatrixExp6(S_mat, thetas(i)) * T1;
    T_configs{i} = T_step;
    plotFrame(T_step, labels{i});
end

%% 3a. Screw axis to [S] Matrix
% https://ethz.ch/content/dam/ethz/special-interest/mavt/robotics-n-intelligent-systems/multiscaleroboticlab-dam/documents/trm/HS2018/Exercise%20Slides/03_2018-10-15_ScrewTheory.pdf
% W5-1 slide 6 
function Smat = screw2Skew(q, s_hat, h)
s_hat = s_hat / norm(s_hat); % norm

% add conditions for pure translation / rotation here? 
% v = -s x q + h * s
v  = cross(-s_hat, q) + h * s_hat;

% se3 
Oskew = [0, s_hat(3), -s_hat(2); -s_hat(3), 0, s_hat(1);s_hat(2), -s_hat(1), 0 ];
Smat = [Oskew, v(:); 0 0 0 0];
end 

%% 3b. Matrix Exp for SE(3)
% https://www.mathworks.com/matlabcentral/answers/1845308-how-to-do-coordinate-transformation-around-a-fixed-axis-using-robotics-toolbox-or-spatial-math-toolb
function T = MatrixExp6(Smat, theta)
    omega_skew = Smat(1:3, 1:3);
    v = Smat(1:3, 4);
    omg = [omega_skew(3,2); omega_skew(1,3); omega_skew(2,1)]; % extract vector
    
    if norm(omg) < 1e-6 % Pure translation w5-1 slide 10
        T = [eye(3), v * theta; 0 0 0 1];
    else
        R = axisangle2rot(omg, theta); % 2a function this is correct from slide 10 w5 - 1 
        % G(theta) = I*theta + (1-cos(theta))[w] + (theta-sin(theta))[w]^2
        G = eye(3)*theta + (1-cos(theta))*omega_skew + (theta-sin(theta))*(omega_skew^2);
        T = [R, G*v; 0 0 0 1];
    end
end

%% 3 Matrix log w5 - 1 slide 12 

function [Smat, theta] = MatrixLog6(T)
% 2 cases 
% if identity no rotation 
end 

%% helper function to plot styuff 
function plotFrame(T, label)
    origin = T(1:3, 4);
    R = T(1:3, 1:3);
    colors = ['r', 'g', 'b'];
    axis_labels = {'x', 'y', 'z'};
    for i = 1:3
        quiver3(origin(1), origin(2), origin(3), R(1,i), R(2,i), R(3,i), 0.5, colors(i), 'LineWidth', 2);
    end
    text(origin(1), origin(2), origin(3), label);
end