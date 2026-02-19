close all; clear all; clc;

%% 1. Define T_sa and Constants
% Origin of {a} relative to {s}
p_sa = [3; 0; 0]; 
% Rotation: x_a=(0,0,1), y_a=(-1,0,0), z_a=(0,1,0)
R_sa = [0 -1  0;
        0  0  -1;
        1  0  0];
T_sa = [R_sa, p_sa; 0 0 0 1];

%% 2. Calculate Screw Parameters (Part i)
theta = acos((trace(R_sa) - 1) / 2);
w_bracket = (1 / (2 * sin(theta))) * (R_sa - R_sa');
w_vec = [w_bracket(3,2); w_bracket(1,3); w_bracket(2,1)];

G_inv = (1/theta)*eye(3) - 0.5*w_bracket + ...
        ((1/theta) - 0.5*cot(theta/2))*(w_bracket^2);
v = G_inv * p_sa;

% Calculate point q on axis and pitch h
h = 0.827;
q = [-1; 1; 0];

%% 3. Plotting
figure('Color', 'w'); hold on; grid on; axis equal;
view(135, 30); xlabel('X_s'); ylabel('Y_s'); zlabel('Z_s');
title('Screw Axis S and Frame {a} in Fixed Frame {s}');

% --- Plot Fixed Frame {s} ---
quiver3(0,0,0, 1,0,0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Xs
quiver3(0,0,0, 0,1,0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Ys
quiver3(0,0,0, 0,0,1, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5); % Zs
text(0,0,0, ' {s}', 'FontSize', 12, 'FontWeight', 'bold');

% --- Plot Frame {a} ---
% Xa (Red), Ya (Green), Za (Blue)
quiver3(p_sa(1), p_sa(2), p_sa(3), R_sa(1,1), R_sa(2,1), R_sa(3,1), 'r--');
quiver3(p_sa(1), p_sa(2), p_sa(3), R_sa(1,2), R_sa(2,2), R_sa(3,2), 'g--');
quiver3(p_sa(1), p_sa(2), p_sa(3), R_sa(1,3), R_sa(2,3), R_sa(3,3), 'b--');
text(p_sa(1), p_sa(2), p_sa(3), ' {a}', 'FontSize', 12);

% --- Plot Screw Axis S ---
% The axis passes through point q and points in direction w_vec
line_len = 3; % length for visualization
axis_line = [q - line_len*w_vec, q + line_len*w_vec];
plot3(axis_line(1,:), axis_line(2,:), axis_line(3,:), 'k', 'LineWidth', 3);
plot3(q(1), q(2), q(3), 'ko', 'MarkerFaceColor', 'k'); % Point q
text(q(1), q(2), q(3), '  Point q', 'FontSize', 10);

% --- Plot Screw Direction Arrow ---
quiver3(q(1), q(2), q(3), w_vec(1), w_vec(2), w_vec(3), 2, 'Color', [1 0.5 0], 'LineWidth', 4);
legend('X_s', 'Y_s', 'Z_s', 'X_a', 'Y_a', 'Z_a', 'Screw Axis', 'q', 's direction');