%% Bicycle simple second-order closed-loop transfer function
% Transfer function: roll angle phi(s) / steer torque T(s)

% Parameters (replace with your own values!)
m  = 94.5;       % mass [kg]
g  = 9.81;     % gravity [m/s^2]
h  = 0.7;    % CoM height [m]
b  = 1.02;      % wheelbase [m]
a  = 0.9;    % horizontal CoM offset
V  = 10;      % forward speed [m/s]
J  = m*h^2;    % inertia approx
D  = a*m*h;    % inertia product
lambda = pi/10; % head angle (rad)
c =0.08;

% Fork / rider parameters (example values)
%k1 = 0.5;      % steer torque-to-angle gain
%k2 = 1.0;      % lean-to-steer feedback gain
k1= b*b/((V*V*sin(lambda)-b*g*cos(lambda))*m*a*c*sin(lambda));  % steer torque-to-angle gain
k2= b*g/(V*V*sin(lambda)-b*g*cos(lambda)); % lean-to-steer feedback gain

% Build transfer function (Control System Toolbox)
num = k1 * [D*V/b, m*V^2*h/b]; % numerator coefficients
den = [J, (D*V/b)*k2, (m*V^2*h/b)*k2 - m*g*h]; % denominator coefficients

G_phi_T = tf(num, den);

disp('Transfer Function G_phi_T(s) = Phi(s)/T(s):');
G_phi_T

%% Define input signal (steer torque over time)
t = 0:0.01:5;% time vector [s]
K = (m*g*a*c*sin(lambda));
u = K*ones(size(t));      % step input torque = 0.1 Nm
% Example: try u = sin(2*pi*1*t) for sinusoidal input

%% Simulate output using lsim
[y, t_out, x] = lsim(G_phi_T, u, t);

%% Plot input vs output
figure;
subplot(2,1,1);
plot(t, u, 'b', 'LineWidth', 1.5);
grid on;
ylabel('Steer Torque Input T(t) [Nm]');
title('Bicycle Response: Steer Torque vs Roll Angle');

subplot(2,1,2);
plot(t_out, y, 'r', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]');
ylabel('Roll Angle \phi(t) [rad]');
