%% THA 1 - PA 
clear all; close all; clc;
% need comments !!! and equations !! add question numbers 

%% rotational matrix SO3 -> equivalent axis angle representation
function [angle, axis] = rot2axisangle(R)
% check if so3
if isSO3(R) == 0
    error("Not a valid SO3 matrix");
end

% calc angle 
cosTheta = (trace(R) - 1) / 2;
angle = acos(cosTheta);
% singularity conditions and then axis calc
if angle < 1e-6
    error("Axis is undefined");
elseif abs(angle - pi) < 1e-6
    sqrtpart = 1/(sqrt(2*(1+R(1,1))));
    axis = [(sqrtpart*(R(1,1)+1)); (sqrtpart*(R(2,1))); (sqrtpart*(R(3,1)))];
    % Correcting signs for the pi case CHECK THIS 
    % if R(1,2) < 0, axis(2) = -axis(2); end
    % if R(1,3) < 0, axis(3) = -axis(3); end
else 
    n = (2 * sin(angle));
    n2 = [R(3,2)-R(2,3), R(1,3) - R(3,1),R(2,1) - R(1,2)]; % do this the easier way maybe?
    axis = n2 / n;
end
end

%% rotational matrix -> Quaternion

function [q0, q1, q2, q3] = rot2quat(R)
% check if so3
if isSO3(R) == 0
    error("Not a valid SO3 matrix");
end
% calculate quat
tr = trace(R);
if (tr>0)
    S = sqrt(tr + 1) * 2; % bro look at this here check sign function again
    % check all the if else statements
    % please write comments here 
    q0 = 0.25 * S; 
    q1 = (R(3,2)-R(2,3)) /S;
    q2 = (R(1,3)-R(3,1)) /S;
    q3 = (R(2,1)-R(1,2)) /S;
end
end
 

%% rotation matrix -> ZYZ

function [phi, theta, psi] = rot2zyz(R)
if isSO3(R) == 0
    error("Not a valid SO3 matrix");
end 
% add singularity
theta = atan2(sqrt(R(1,3)^2 + R(2,3)^2), R(3,3));
% need to add inbetween -pi and 0 whoops; fix statements
if abs(theta) < 1e-6 || abs(theta - pi) < 1e-6
    phi = 0;
    if R(3,3) > 0
        theta = 0; % infinity solutions throw error for pi, -pi and 0 
        psi = atan2(-R(1,2), R(1,1));
    else 
        theta = pi;
        psi = atan2(R(1,2), -R(1,1));
    end 
else 
    phi = atan2(R(2,3),R(1,3));
    psi = atan2(R(3,2),-R(3,1));
end 
end  


%% rotation matrix -> roll pitch yaw

function [roll, pitch, yaw] = rot2zyx(R)
if isSO3(R) == 0
    error("Not a valid SO3 matrix");
end 
% add singularity and change to x y z
pitch = atan2(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2));

if abs(pitch - pi/2) < 1e-6 || abs(pitch + pi/2) < 1e-6
    yaw = 0;
    if pitch > 0
        roll = atan2(R(1,2), R(2,2));
    else 
        roll = -atan2(R(1,2), R(2,2));
    end 
else 
    roll = atan2(R(3,2), R(3,3));
    yaw = atan2(R(2,1), R(1,1));
end 
end 
%% axis angle -> rotation

function [R] = axisangle2rot(axis, angle)

% unit vector check
if norm(axis)>0
    axis = axis / norm(axis);
end 

ux = axis(1);
disp(ux);
uy = axis(2);
disp(uy);
uz = axis(3);
disp(uz);

K = [0, -uz, uy; uz, 0, -ux; -uy, ux, 0];
disp("K")
disp(K);

% rodriguez formula 
R = eye(3) + sin(angle)*K + (1 - cos(angle))*K^2; % I + wsintheta + w2 (1-costheta)

end 

%% quat -> rotation
function [R] = quat2rot(q0, q1, q2, q3)
% ALL THE SIGNS R WRONGGGGGG 
% normalized here dk why 
mag = sqrt(q0^2 + q1^2 + q2^2 + q3^2);
q0 = q0/mag; q1 = q1/mag; q2 = q2/mag; q3 = q3/mag;

r11 = 1 - 2*q2^2 - 2*q3^2;
r12 = 2*q1*q2 - 2*q0*q3;
r13 = 2*q1*q3 + 2*q0*q2;

r21 = 2*q1*q2 + 2*q0*q3;
r22 = 1 - 2*q1^2 - 2*q3^2;
r23 = 2*q2*q3 - 2*q0*q1;

r31 = 2*q1*q3 - 2*q0*q2;
r32 = 2*q2*q3 + 2*q0*q1;
r33 = 1 - 2*q1^2 - 2*q2^2;

R = [r11, r12, r13;
    r21, r22, r23;
    r31, r32, r33];
end


%% Rotation matrix check

function SO3 = isSO3(R)

% 3 by 3 matrix
is3by3 = all(size(R) == [3 3]);

% orthogonal: R'R = I
isOrthog = norm(R' * R - eye(3)) < 1e-6;

%  determinant: det(R) = 1
isDetOne = abs(det(R) - 1) < 1e-6;

SO3 = is3by3 && isOrthog && isDetOne;
end

%% Main / Testing functions

R = [0 -1 0; 1 0 0; 0 0 1]; % 90 degree rotation around Z

[angle, axis] = rot2axisangle(R);
[R2] = axisangle2rot(axis, angle);
[q0, q1, q2, q3] = rot2quat(R);
[R3] = quat2rot(q0, q1, q2, q3);
[phi, theta, psi] = rot2zyz(R);
[roll, pitch, yaw] = rot2zyx(R);

disp("Rotational Matrix:"); disp(R);
fprintf("Angle: %.4f rad\n", angle);
fprintf("Axis: [%.4f, %.4f, %.4f]\n", axis);
fprintf("Quaternion: [%.4f, %.4f, %.4f, %.4f]\n", q0, q1, q2, q3);
fprintf("ZYZ: [%.4f, %.4f, %.4f]\n", phi, theta, psi);
fprintf("ZYX: [%.4f, %.4f, %.4f]\n", roll, pitch, yaw);
disp("/n");

disp("From Axis Angle")
disp("Rotational Matrix:"); disp(R2);
disp("From Quaternion")
disp("Rotational Matrix:"); disp(R3);

% Checking with built in matlab functions :P 
disp("Checks with builtin matlab functions");
disp("Rotation Matrix 2 Axis Angle"); disp(rotm2axang(R));
disp("Rotation Matrix 2 Quat"); disp(rotm2quat(R));
disp("Rotation Matrix 2 ZYZ"); disp(rotm2eul(R, 'ZYZ'));
disp("Rotation Matrix 2 ZYX"); disp(rotm2eul(R, 'XYZ'));



