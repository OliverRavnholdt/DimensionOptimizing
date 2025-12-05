function [Phi_system, Phi_q, nu, gamma] = setup_constraints()
%% Definitions
% Define angles
beta = 5/4 * pi;

% Define lengths
syms L2 L32 L4 L51 L52 L6 alpha_sym
L1 = 24e-3;
L31 = 2/3 * L32 *sin(alpha_sym)/alpha_sym; % m Centroid without holes and using joint lengths
L53 = L51/2; % m Centroid without holes and usin joint length
SL = 22e-3; % m
SH = 20e-3; % m
L = dictionary(["L2", "L32", "L4", "L51", "L52", "L6", "alpha"], [L2 L32, L4, L51, L52, L6, alpha_sym]);

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
q = lib2D.getq(x,y,phi,nb);
q_dot = lib2D.getq(x_dot,y_dot,phi_dot,nb);

% Define orientation matrices
A = cell(1, nb);
for i = 1:nb
    A{i} = lib2D.orientationMatrix2D(phi(i));
end

% Define r vectors
r = cell(1, nb);
for i = 1:nb
    r(i) = {[x(i); y(i)]};
end

% Define S'^Q vectors
S1_prime_Q = [-L1/2; 0];
S3_prime_Q = lib2D.orientationMatrix2D(beta)*[L31; 0];
S6_prime_Q = [-L6/2; 0];

% Define S'^P vectors
S1_prime_P1 = [L1/2; 0];
S2_prime_P2 = [-L2/2; 0];
S2_prime_P3 = [L2/2; 0];
S3_prime_P4 = S3_prime_Q+lib2D.orientationMatrix2D(alpha_sym-pi/4)*[0; L32];
S3_prime_P5 = S3_prime_Q+lib2D.orientationMatrix2D(pi/4-alpha_sym)*[L32; 0];
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

Phi_D(1) = phi(1) - omega1*t+pi;
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
end