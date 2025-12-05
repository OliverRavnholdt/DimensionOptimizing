% Function that returns all omega pairs for given grid size without 0
% rotation pairs
function omega_pairs = get_omega_pairs(grid_size)
omega1_set = round(linspace(-1, 1, grid_size) * pi,2);
omega6_set = round(linspace(-1, 1, grid_size) * pi,2);

omega_pairs = combvec(omega1_set, omega6_set)';  % all combinations

% Remove 0 rotation entries
omega_pairs = omega_pairs(~(omega_pairs(:,1) == 0 & omega_pairs(:,2) == 0), :);
end