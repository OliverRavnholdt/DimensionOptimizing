% Keep codespace clean to avoid bugs and improve readability
clear; clc;

%% Definitions
% Define angles
beta = 5/4 * pi;

% Define lengths
syms L2 L32 L4 L51 L52 L6 alpha
L1 = 24e-3;
L31 = 2/3 * L32 *sin(alpha)/alpha; % m Centroid without holes and using joint lengths
L53 = L51/2; % m Centroid without holes and using joint length
SL = 22e-3; % m
SH = 20e-3; % m
L = dictionary(["L2", "L32", "L4", "L51", "L52", "L6", "alpha"], [L2 L32, L4, L51, L52, L6, alpha]);

% Define number of bodies
nb = 6;

% Define unknowns
syms x [1 nb]
syms y [1 nb]
syms phi [1 nb]
syms x_dot [1 nb]
syms y_dot [1 nb]
syms phi_dot [1 nb]
syms t

% Define q and q_dot
q = lib.getq(x,y,phi,nb);
q_dot = lib.getq(x_dot,y_dot,phi_dot,nb);
%%
% Define orientation matrices
A = cell(1, nb);
for i = 1:nb
    A{i} = lib.orientationMatrix2D(phi(i));
end

% Define r vectors
r = cell(1, nb);
for i = 1:nb
    r(i) = {[x(i); y(i)]};
end

% Define S'^Q vectors
S1_prime_Q = [-L1/2; 0];
S3_prime_Q = lib.orientationMatrix2D(beta)*[L31; 0];
S6_prime_Q = [-L6/2; 0];

% Define S'^P vectors
S1_prime_P1 = [L1/2; 0];
S2_prime_P2 = [-L2/2; 0];
S2_prime_P3 = [L2/2; 0];
S3_prime_P4 = S3_prime_Q+lib.orientationMatrix2D(alpha-pi/4)*[0; L32];
S3_prime_P5 = S3_prime_Q+lib.orientationMatrix2D(pi/4-alpha)*[L32; 0];
S4_prime_P6 = [-L4/2; 0];
S4_prime_P7 = [L4/2; 0];
S5_prime_P8 = [-L53; 0];
S5_prime_P9 = [-(L53-L52); 0];
S6_prime_P10 = [L6/2; 0];

% Define C vectors
C1 = [0; 0];
C2 = [SL; SH];

%% Constraints
% Define kinematic constraints
Phi_K = sym(zeros(16, 1));

% Absolute Constraints
Phi_K(1:2) = r(1) + A{1}*S1_prime_Q - C1;
Phi_K(3:4) = r(3) + A{3}*S3_prime_Q - C2;
Phi_K(5:6) = r(6) + A{6}*S6_prime_Q - C2;

% Relative Constraints
Phi_K(7:8) = r(2) + A{2}*S2_prime_P2 - r(1) - A{1}*S1_prime_P1;
Phi_K(9:10) = r(3) + A{3}*S3_prime_P4 - r(2) - A{2}*S2_prime_P3;
Phi_K(11:12) = r(4) + A{4}*S4_prime_P6 - r(3) - A{3}*S3_prime_P5;
Phi_K(13:14) = r(5) + A{5}*S5_prime_P8 - r(4) - A{4}*S4_prime_P7;
Phi_K(15:16) = r(6) + A{6}*S6_prime_P10 - r(5) - A{5}*S5_prime_P9;

% Define driving constraint
Phi_D = sym(zeros(2,1));

syms omega1 omega6

Phi_D(1) = phi(1) - omega1*t-pi+pi/4;
Phi_D(2) = phi(6) - omega6*t+1/4*pi;

% Define Phi system
Phi_system = [Phi_K; Phi_D];

%% Get jacobian, nu, gamma in matlabFunctions
% Define jacobian using MatLab built in function
Phi_q = jacobian(Phi_system, q);

% Define nu
nu = -diff(Phi_system, t);

% Define gamma
gamma = -jacobian(Phi_q*q_dot, q)*q_dot-2*diff(Phi_q, t)*q_dot+diff(nu,t);

% Convert from symbolic functions to matlabFunctions
x_sym = [values(L); omega1; omega6];
Phi_system = matlabFunction(Phi_system, 'Vars', {q, t, x_sym});
nu = matlabFunction(nu, 'Vars', {q, x_sym});
gamma = matlabFunction(gamma, 'Vars', {q, q_dot, t, x_sym});
Phi_q = matlabFunction(Phi_q, 'Vars', {q, x_sym});


%% Rotation Grid
omega1_set = round(linspace(-1, 1,5) * pi,2);
omega6_set = round(linspace(-1, 1, 5) * pi,2);

omega_pairs = combvec(omega1_set, omega6_set)';  % all combinations
omega_pairs = omega_pairs(~(omega_pairs(:,1) == 0 & omega_pairs(:,2) == 0), :);

%% NR setup
Q_all = cell(size(omega_pairs,1), 1);

% Options
options = struct();
options.normCorrection = 1.0;
options.tolerance = 1e-12;
options.maxIter = 1000;

% Initial guess of lengths that works
alpha_0 = pi/4;
L2_0 = 25e-3; % m
L32_0 = 34.5e-3; % m
L31_0 = 2/3 * L32_0 *sin(alpha_0)/alpha_0; % m Centroid without holes and using joint lengths
L4_0 = 97e-3; % m
L51_0 = 111e-3; % m
L52_0 = 31.5e-3; % m
L53_0 = L51_0/2; % m Centroid without holes and usin joint length
L6_0 = 89e-3; % m
L0 = [L2_0, L32_0, L4_0, L51_0, L52_0, L6_0, alpha_0]';

% Newton raphson solver for each omega pair
tic % Time start
for i = 1:size(omega_pairs,1)
% Define current omega pairs
omega0 = [omega_pairs(i, 1), omega_pairs(i,2)]';

% Define column vector of initial guess and omega pair
x0 = [L0; omega0];

% Angles meassured from sketch in Solidworks
angles = [pi-pi/4, 0.96, 0.925, 2*pi-0.65, pi+0.5, 2*pi-0.61+1];

% Get position of all bodies as initial guess
Q0 = subs([L1/2*cos(angles(1)), L1/2*sin(angles(1)), angles(1), ...
          L1*cos(angles(1))+L2/2*cos(angles(2)), L1*sin(angles(1))+L2/2*sin(angles(2)), angles(2), ...
          SL+L31*cos(angles(3)), SH+L31*sin(angles(3)), angles(3), ...
          SL+cos(angles(6))*L6+L52*cos(angles(5))+L4/2*cos(angles(4)), SH+sin(angles(6))*L6+L52*sin(angles(5))+L4/2*sin(angles(4)), angles(4), ...
          SL+cos(angles(6))*L6+(L53-L52)*cos(angles(5)), SH+sin(angles(6))*L6+(L53-L52)*sin(angles(5)), angles(5), ...
          SL+cos(angles(6))*L6/2, SH+sin(angles(6))*L6/2, angles(6)], values(L), L0);

% Matrix of all time steps
d0 = 0;
dt = 100;
% Try to rotate largest omega by 180 deg
t_end = pi / max(abs(omega0));
time = linspace(d0, t_end, dt);
%% NR solver
% Define the wanted lengths
Lw = [0.029964840164286	0.027764051138324	0.101459898079935	0.131616339441701	0.037731970159884	0.080320029247764	0.718887343640291]';

% Get difference from initial guess to wanted lengths
dL = Lw-L0;

% Prepare gradual length changes
nSteps = 20; % Steps

% Define matrix of length path to solve to get final wanted length
L_path = L0 + dL .* linspace(0, 1, nSteps+1);

% Define column vector of path to final length and omega pair
x_path = [L_path; omega0*ones(1,size(L_path,2))]; % Define

%% Change lengths gradually
% Exact initial solution from initial guess using NR soliver
[feasable, Q] = lib.NRLengthChange(Q0, Phi_system, Phi_q, nu, gamma, time, nb, x0, options);
Q0 = Q(:,1);

% Solve the start position of each length in path
for j = 2:size(L_path, 2)  % start from second step (since first = initial)
    x_curr = x_path(:, j);  % current lengths

    % Run Newton–Raphson solver
    [feasable, Q] = lib.NRLengthChange(Q0, Phi_system, Phi_q, nu, gamma, 0, nb, x_curr, options);

    % If solver failed, you can add diagnostics
    if ~feasable
        fprintf('NR failed at step %d/%d\n', j, size(x_curr, 2));
        break
    end

    % Update for next iteration
    Q0 = Q(:, 1);
end

% Solve the entire system for all time steps for wanted length
[feasable, Q] = lib.NRLengthChange(Q0, Phi_system, Phi_q, nu, gamma, time, nb, x_curr(:, end), options);

% Make sure that solver did not find angle that was flipped as this could
% be a possible solution in some cases. If it did the following code will
% compare the angle with the previous angle and map the 180deg rotation to
% the closest new angle to make sure unwanted flips does not happen
nT = size(Q,2);
for k = 2:nT
    for bi = 1:nb
        % Get angles
        phi_idx = 3*bi;
        phi_prev = Q(phi_idx, k-1);
        phi_curr = Q(phi_idx, k);
        % Map the current angle to the nearest equivalent of phi_prev
        d = wrapToPi(phi_curr - phi_prev);   % requires Mapping Toolbox
        Q(phi_idx, k) = phi_prev + d;
    end
end
Q_all{i} = Q;
end
toc % Time stop

%% Animate movement
% Setup
fig1 = figure; hold on; axis square;
xlim([-0.2,0.2])
ylim([-0.2,0.2])
title("Leg Animation for all omega pairs")
xlabel("X (mm)")
ylabel("Y (mm)")

colors = ["#fcba03", "#b3091d", "#347a05", "#05637a", "#670b9c", "#a6144e"];
legends = ["Body1", "Body2", "Body3", "Body4", "Body5", "Body6"];

b1plot = plot(nan(1,2), nan(1,2), "Color", colors(1));
b2plot = plot(nan(1,2), nan(1,2), "Color", colors(2));
b3plot1 = plot(nan(1,2), nan(1,2), "Color", colors(3));
b3plot2 = plot(nan(1,2), nan(1,2), "Color", colors(3), "HandleVisibility","off");
b3plot3 = plot(nan(1,2), nan(1,2), "Color", colors(3), "HandleVisibility","off");
b3plot4 = plot(nan(1,2), nan(1,2), "Color", colors(3), "HandleVisibility","off");
b3plot5 = plot(nan(1,2), nan(1,2), "Color", colors(3), "HandleVisibility","off");
b4plot = plot(nan(1,2), nan(1,2), "Color", colors(4));
b5plot = plot(nan(1,2), nan(1,2), "Color", colors(5));
b6plot = plot(nan(1,2), nan(1,2), "Color", colors(6));
legend(legends)

% Loop over all omega pairs
for j = 1:size(omega_pairs, 1)
% Extract the corresponding Q
Q = Q_all{j};

% Define S'^Q vectors
L0 = L_path(:,end);
S1_prime_Q = [-L1/2; 0];
S3_prime_Q = lib.orientationMatrix2D(beta) * [2/3 * L0(2) * sin(L0(7)) / L0(7); 0];
S6_prime_Q = [-L0(6)/2; 0];

% Define S'^P vectors
S1_prime_P1 = [L1/2; 0];
S2_prime_P2 = [-L0(1)/2; 0];
S2_prime_P3 = [L0(1)/2; 0];
S3_prime_P4 = S3_prime_Q + lib.orientationMatrix2D(L0(7)-pi/4)*[0; L0(2)];
S3_prime_P5 = S3_prime_Q + lib.orientationMatrix2D(pi/4-L0(7))*[L0(2); 0];
S4_prime_P6 = [-L0(3)/2; 0];
S4_prime_P7 = [L0(3)/2; 0];
S5_prime_P8 = [-L0(4)/2; 0];
S5_prime_P9 = [-(L0(4)/2 - L0(5)); 0];
S6_prime_P10 = [L0(6)/2; 0];

% Define body geometry from COM
body1 = [Q(1,:) - L1/2 * cos(Q(3,:));  Q(2,:) - L1/2 * sin(Q(3,:));
         Q(1,:) + L1/2 * cos(Q(3,:));  Q(2,:) + L1/2 * sin(Q(3,:))];

body2 = [Q(4,:) - L0(1)/2 * cos(Q(6,:));  Q(5,:) - L0(1)/2 * sin(Q(6,:));
         Q(4,:) + L0(1)/2 * cos(Q(6,:));  Q(5,:) + L0(1)/2 * sin(Q(6,:))];

body3abs = zeros(2, length(time));
body3P5  = zeros(2, length(time));
body3P4  = zeros(2, length(time));
P4test   = zeros(2, length(time));
P5test   = zeros(2, length(time));

for i = 1:size(Q,2)
    body3abs(:,i) = Q([7,8],i) + lib.orientationMatrix2D(Q(9,i)) * S3_prime_Q;
    body3P5(:,i)  = Q([7,8],i) + lib.orientationMatrix2D(Q(9,i)) * S3_prime_P5;
    body3P4(:,i)  = Q([7,8],i) + lib.orientationMatrix2D(Q(9,i)) * S3_prime_P4;
end

body4 = [Q(10,:) - L0(3)/2 * cos(Q(12,:));  Q(11,:) - L0(3)/2 * sin(Q(12,:));
         Q(10,:) + L0(3)/2 * cos(Q(12,:));  Q(11,:) + L0(3)/2 * sin(Q(12,:))];

body5 = [Q(13,:) - L0(4)/2 * cos(Q(15,:));  Q(14,:) - L0(4)/2 * sin(Q(15,:));
         Q(13,:) + L0(4)/2 * cos(Q(15,:));  Q(14,:) + L0(4)/2 * sin(Q(15,:))];

body6 = [Q(16,:) - L0(6)/2 * cos(Q(18,:));  Q(17,:) - L0(6)/2 * sin(Q(18,:));
         Q(16,:) + L0(6)/2 * cos(Q(18,:));  Q(17,:) + L0(6)/2 * sin(Q(18,:))];

% Plot C1 and C2
plot(C1(1), C1(2), "ko", "HandleVisibility","off")
plot(C2(1), C2(2), "ko", "HandleVisibility","off")

% Let animation render
pause(1)

% Animate plot
for i = 1:size(Q,2)
    b1plot.XData = [body1([1,3], i)];
    b1plot.YData = [body1([2,4], i)];

    b2plot.XData = [body2([1,3], i)];
    b2plot.YData = [body2([2,4], i)];

    b3plot1.XData = [body3P4(1,i), body3P5(1,i)];
    b3plot1.YData = [body3P4(2,i), body3P5(2,i)];

    b3plot2.XData = [body3abs(1,i), body3P5(1,i)];
    b3plot2.YData = [body3abs(2,i), body3P5(2,i)];

    b3plot3.XData = [body3abs(1,i), body3P4(1,i)];
    b3plot3.YData = [body3abs(2,i), body3P4(2,i)];

    b4plot.XData = [body4([1,3], i)];
    b4plot.YData = [body4([2,4], i)];

    b5plot.XData = [body5([1,3], i)];
    b5plot.YData = [body5([2,4], i)];

    b6plot.XData = [body6([1,3], i)];
    b6plot.YData = [body6([2,4], i)];

    drawnow;
    pause(0.05)
end
end


%% Plot kinematics
Q = Q_all{end};
Q_dot = zeros(nb*3, size(Q,2));
Q_ddot = zeros(nb*3, size(Q,2));
for i = 1:size(Q,2)
    Q_dot(:,i) = Phi_q(Q(:,i), x0) \ nu(Q(:,i), x0);
    Q_ddot(:,i) = Phi_q(Q(:,i), x0) \ gamma(Q(:,i), Q_dot(:,i), time(i), x0);
end

time_plot = time(:, 1:size(Q,2));
fig2 = lib.plotKinematics(Q,Q_dot,Q_ddot, time_plot, nb);
