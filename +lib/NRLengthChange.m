% Newton raphson solver that returns Q: A matrix of all positions and
% angles for each timestep and returns if the solver found a physically
% viable solution based on custom tuned parameters and the convergence of
% the solver itself
function [converged, Q] = NRLengthChange(q0, Phi_system, Phi_q, nu, gamma, time, nb, x0, options)
    arguments
       q0;
       Phi_system;
       Phi_q;
       nu,
       gamma,
       time;
       nb double;
       x0;
       options struct = struct()
    end

% Setup
Q = zeros(nb*3, length(time));
Q_dot = zeros(nb*3, size(Q,2));
Q_ddot = zeros(nb*3, size(Q,2));
angular_acc = zeros(nb, size(Q,2));
converged = true;

% Initial position
Q(:,1) = q0;

% Loop Newton Raphson for all timesteps
for i = 1:length(time)
    % Make initial guess the found value from last timestep
    if i ~= 1
        Q(:,i) = Q(:,i-1);
    end

    % start improving the guess
    normCorrection = options.normCorrection;
    tolerance = options.tolerance;
    iterationCounter = 1;
    maxIter = options.maxIter;
    while normCorrection>tolerance && iterationCounter <= maxIter
        residual = Phi_system(Q(:,i), time(i), x0);
        jacobian = Phi_q(Q(:,i), x0);
        correction = jacobian \ residual;
        Q(:,i) = Q(:,i) - correction;
        normCorrection = norm(correction);
        iterationCounter = iterationCounter + 1;
    end

% Get accelerations and velocities
Phi_q_curr = Phi_q(Q(:,i), x0);
Q_dot(:,i) = Phi_q_curr \ nu(Q(:,i), x0);
Q_ddot(:,i) = Phi_q_curr \ gamma(Q(:,i), Q_dot(:,i), time(i), x0);

angular_acc(:,i) = [Q_ddot(3:3:end,i)];
max_angular_acc = max(abs(angular_acc(:,i)));

% Different parameters that have been manually tuned to make sure that
% the system does not rotate bodies more than physically posssible

% Angle of 1 and 2 cannot be equal and 2<1
cond_angle12 = max_angular_acc > 100 || mod(Q(3,i), 2*pi) - mod(Q(6,i), 2*pi) <= 0;

% Angle 3 and 4 should not be rotated 180 deg
cond_angle34 = abs(Q(9,i) - (mod(Q(12,i), 2*pi)-pi)) < 0.3;

% Angle 5 and 6 should not be rotated 180 deg
cond_angle56 = abs((mod(Q(15,i), 2*pi)-pi) - mod(Q(18,i), 2*pi)) < 0.3;

% Angle 1 should have upper limit
cond_angle1 = mod(Q(3,i), 2*pi) > pi+pi/4;

% Angle 6 and 3 should not be rotated 180 deg
cond_angle36 = abs(mod(Q(9,i), 2*pi) - mod(Q(18,i), 2*pi)) < 0.3;

% NR should converge
Cond_converge = normCorrection > tolerance;

% All conditions should be false
cond_eval = cond_angle1 + cond_angle12 + cond_angle34 + cond_angle56 + cond_angle36 + Cond_converge;

    % Check if system is good
    if  cond_eval ~= 0
        % warning("Newton-Raphson did not converge at timestep %d", i);
        converged = false;
        Q = Q(:,1:i);
        return;
    end
end
end