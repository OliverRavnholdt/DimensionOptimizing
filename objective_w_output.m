% Cost function with extra data output for comparing data from ga solution
function out = object_w_output(Q0, Phi_system, Phi_q, nu, gamma, t, nb, omegas, Lw, L0, wpos, options)
out = struct();

% Copy inputs as this script is setup to potentially run parfar loop
L_cur = L0(:);
L_target = Lw(:);

% Setup
SL = 20e-3; % m
Q_all = cell(1,length(omegas));

% Loop over all omega pairs for better explanation see
% KinematicsLengthChange
for i = 1:length(omegas)
% Iniduvidual setup for each omega pair
Q = [];
x_curr = [];
Q0_local = Q0;
omega_row = omegas(i, :);
omega0 = omega_row(:);
t_local = t;
t_end = pi / max(abs(omega_row));
t_local.end = t_end;
min_foot_height = 0;
max_foot_height = -20e6;

time = linspace(t_local.t0, t_local.end, t_local.steps);

%% Prepare gradual length changes
dL = L_target-L_cur;

nSteps = 20;
L_path = L_cur + dL .* linspace(0, 1, nSteps+1);
x_path = [L_path; omega0*ones(1,size(L_path,2))];

for j = 2:size(L_path, 2)
    x_curr = x_path(:, j);

    % Run Newton–Raphson solver
    [feasable, Q] = lib.NRLengthChange(Q0_local, Phi_system, Phi_q, nu, gamma, 0, nb, x_curr, options);

    % If solver failed punish cost as wanted length cannot be reached 
    if ~feasable
        % Huge cost to not use this solution
        cost = 1e6;
        return
    end

    % Update for next iteration
    Q0_local = Q(:, 1);
end

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

% Final solution for all timesteps
[~, Q] = lib.NRLengthChange(Q0_local, Phi_system, Phi_q, nu, gamma, time, nb, x_curr(:, end), options);

% Get the foot position for all timesteps
try
Foot_pos = [Q(13,:) + L0(4)/2 * cos(Q(15,:));  Q(14,:) + L0(4)/2 * sin(Q(15,:))];

% Find the tallest and lowest point of foot position in y.direction while
% foot is between servo motors in x-direction
for k = 1:length(Foot_pos)
    curr_foot_pos = Foot_pos(:,k);
if curr_foot_pos(1) > 0 && curr_foot_pos(1) < SL
    min_foot_height = min(curr_foot_pos(2), min_foot_height);
    if curr_foot_pos(2) <= 0
        max_foot_height = max(curr_foot_pos(2), max_foot_height);
    end
end
end
catch
end

% Store all Q values for all omega pairs
Q_all{i} = Q;

end
% Cost weights
w_total_pos = 1;
w_pos_upper = 1e4;
w_pos_lower = 1e4;

% Cost
% First the cost is the total timesteps for all omega pairs
cost = t.steps * length(omegas);
% Positions where the NR could converge that are physcially possible
total_pos = sum(cellfun(@(x) size(x, 2), Q_all));

% Subtract positions that converged from cost as well as the lowest point.
% Add the tallest point to cost.
cost = cost - w_total_pos * total_pos + w_pos_upper*abs(max_foot_height) - w_pos_lower*abs(min_foot_height);

% Extra output data
out.total_pos = total_pos;
out.cost = cost;
out.wpos = wpos;
out.max_foot_height = max_foot_height;
out.min_foot_height = min_foot_height;
end