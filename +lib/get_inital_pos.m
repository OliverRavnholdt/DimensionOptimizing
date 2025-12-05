% Setup function for GA that returns initial position solved with NR and
% initial lengths
function [Q0, L0] = get_inital_pos(Phi_system, Phi_q, nu, gamma, nb, options)
% Define lengths
syms L2 L32 L4 L51 L52 L6 alpha_sym
L1 = 24e-3;
L31 = 2/3 * L32 *sin(alpha_sym)/alpha_sym; % m Centroid without holes and using joint lengths
L53 = L51/2; % m Centroid without holes and usin joint length
SL = 22e-3; % m
SH = 20e-3; % m
L = dictionary(["L2", "L32", "L4", "L51", "L52", "L6", "alpha"], [L2 L32, L4, L51, L52, L6, alpha_sym]);

% Initial guess that works
alpha_0 = pi/4;
L2_0 = 25e-3; % m
L32_0 = 34.5e-3; % m
L4_0 = 97e-3; % m
L51_0 = 111e-3; % m
L52_0 = 31.5e-3; % m
L6_0 = 89e-3; % m
L0 = [L2_0, L32_0, L4_0, L51_0, L52_0, L6_0, alpha_0]';

% Angles meassured from sketch in Solidworks
angles = [pi-pi/4, 0.96, 0.925, 2*pi-0.65, pi+0.5, 2*pi-0.61+1];

% Get position of all bodies as initial guess
Q0 = subs([L1/2*cos(angles(1)), L1/2*sin(angles(1)), angles(1), ...
          L1*cos(angles(1))+L2/2*cos(angles(2)), L1*sin(angles(1))+L2/2*sin(angles(2)), angles(2), ...
          SL+L31*cos(angles(3)), SH+L31*sin(angles(3)), angles(3), ...
          SL+cos(angles(6))*L6+L52*cos(angles(5))+L4/2*cos(angles(4)), SH+sin(angles(6))*L6+L52*sin(angles(5))+L4/2*sin(angles(4)), angles(4), ...
          SL+cos(angles(6))*L6+(L53-L52)*cos(angles(5)), SH+sin(angles(6))*L6+(L53-L52)*sin(angles(5)), angles(5), ...
          SL+cos(angles(6))*L6/2, SH+sin(angles(6))*L6/2, angles(6)], values(L), L0);

% Pass in any rotation as the position is only solved for t=0
x0 = [L0; 0; 0];

% Get initial position using NR solver
[~, Q0] = lib.NRLengthChange(Q0, Phi_system, Phi_q, nu, gamma, 0, nb, x0, options);
end