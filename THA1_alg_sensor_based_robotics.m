%% THA 1 - PA 
clear all; close all; clc;

%% rotational matrix SO3 -> equivalent axis angle representation
function [angle, axis] = rot2axis(R)
% check if so3
if isSO3(R) == 0
    disp("Not a valid rotation matrix");
    return
end

% calc angle 
cosTheta = (trace(R) - 1) / 2;
angle = acos(cosTheta);
% singularity conditions and then axis calc
if angle == 0
    axis = NaN;
elseif angle == pi
    n = 1 / (sqrt(2*(1+R(3,3))));
    n2 = [R(1,3),R(2,3),1+R(3,3)];
    axis = n2 / n;
else 
    n = (1/2 * sin(angle));
    n2 = [R(3,2)-R(2,3), R(1,3) - R(3,1),R(2,1) - R(1,2)];
    axis = n2 / n;
end 

end

%% rotational matrix -> Quaternion

function [q0, q1, q2, q3] = rot2quat(R)
% check if so3
if isSO3(R) == 0
    disp("Not a valid rotation matrix");
    return
end
% calculate quat
q0 = 1/2 * sqrt(R(1,1)+R(2,2)+R(3,3)+1);
q1 = 1/2 * sign((R(3,2)-R(2,3))*sqrt(R(1,1)-R(2,2)-R(3,3)+1));
q2 = 1/2 * sign((R(1,3)-R(3,1))*sqrt(R(2,2)-R(3,3)-R(1,1)+1));
q3 = 1/2 * sign((R(2,1)-R(1,2))*sqrt(R(3,3)-R(1,1)-R(2,2)+1));
end 

%% rotation matrix -> ZYZ



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
[angle, axis] = rot2axis(R);
[q0, q1, q2, q3] = rot2quat(R);

% Displaying rotation 
disp("Rotational Matrix");
disp(R)

% Dispaly axis angle 
disp("Axis");
disp(axis);
disp("Angle");
disp(angle);

% Display Quat
disp("Quat");
disp(q0);
disp(q1);
disp(q2);
disp(q3);


if isSO3(R)
    disp("so3 check is true");
end 

