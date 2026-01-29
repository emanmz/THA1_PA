%% THA 1 - PA 
clear all; close all; clc;

%% Given a rotational matrix SO3 -> equivalent axis angle representation
function [angle, axis] = rotation2axis(R)
    % calculate rotation anglle
    % trace(R) = R(1,1) + R(2,2) + R(3,3)
    cosTheta = (trace(R) - 1) / 2;
    angle = acos(cosTheta);
    n = (1/2 * sin(angle));
    n2 = [R(3,2)-R(2,3), R(1,3) - R(3,1),R(2,1) - R(1,2)];
    axis = n * n2;
        
end

R = [0 -1 0; 1 0 0; 0 0 1]; % 90 degree rotation around Z
[angle, axis] = rotation2axis(R);