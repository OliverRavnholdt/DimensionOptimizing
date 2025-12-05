clear; clc;

clear; clc;
% System
[Phi_system, Phi_q, nu, gamma] = lib.setup_constraints();
nb = 6;

% Options
options = struct();
options.normCorrection = 1.0;
options.tolerance = 1e-8;
options.maxIter = 1000;

% Initial Guess
[Q0, L0] = lib.get_inital_pos(Phi_system, Phi_q, nu, gamma, nb, options);

% Omega grid
grid_size = 5;
omegas = lib.get_omega_pairs(grid_size);

% Time
time = struct();
time.steps = 100;
time.t0 = 0;

wpos = 1e3;
Lw = L0;
out = objective_w_output(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, wpos, options);
names = {'L2','L32','L4','L51','L52','L6','alpha'};
tableData = L0;
tableNames = ["L0"];
fprintf("For inital guess weight: %0.5f: \n" + ...
    "Cost: %0.5f\n" + ...
    "Leg Max Height: %0.5f\n" + ...
    "Leg Min Height: %0.5f\n" + ...
    "Total positions: %i\n\n", wpos, out.cost, out.max_foot_height, out.min_foot_height, out.total_pos);
runs = cell(1,5);
%% Mesh 5 Weight 0
wpos = 0;
Lw = [0.025678407485219	0.028202216934747	0.087564818345634	0.129777496908124	0.033334758164504	0.089457299731085	0.641853909462706]';
out = objective_w_output(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, wpos, options);
fprintf("For inital guess weight: %0.5f: \n" + ...
    "Cost: %0.5f\n" + ...
    "Leg Max Height: %0.5f\n" + ...
    "Leg Min Height: %0.5f\n" + ...
    "Total positions: %i\n\n", wpos, out.cost, out.max_foot_height, out.min_foot_height, out.total_pos);
tableData = [tableData, Lw];
tableNames = [tableNames, "M5 W0"];
%% Mesh 5 Weight e3
wpos = 1e3;
Lw = [0.026149594043088	0.027899330640565	0.089482609277939	0.123808530730718	0.036692269756823	0.089391963635236	0.628934926575700]';
out = objective_w_output(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, wpos, options);
fprintf("For inital guess weight: %0.5f: \n" + ...
    "Cost: %0.5f\n" + ...
    "Leg Max Height: %0.5f\n" + ...
    "Leg Min Height: %0.5f\n" + ...
    "Total positions: %i\n\n", wpos, out.cost, out.max_foot_height, out.min_foot_height, out.total_pos);
tableData = [tableData, Lw];
tableNames = [tableNames, "M5 We3"];
%% Mesh 5 Weight e4
wpos = 1e4;
Lw = [0.027747141213146	0.028223795675887	0.092341309186911	0.133199998285012	0.031261306821719	0.079961223873390	0.880158011192123]';
out = objective_w_output(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, wpos, options);
fprintf("For inital guess weight: %0.5f: \n" + ...
    "Cost: %0.5f\n" + ...
    "Leg Max Height: %0.5f\n" + ...
    "Leg Min Height: %0.5f\n" + ...
    "Total positions: %i\n\n", wpos, out.cost, out.max_foot_height, out.min_foot_height, out.total_pos);
tableData = [tableData, Lw];
tableNames = [tableNames, "M5 We4"];
%% Mesh 5 Weight e5
wpos = 1e5;
Lw = [0.028157061424391	0.031899915329519	0.108126539186208	0.133199999471170	0.025200003691616	0.092255339389248	0.843490428705472]';
out = objective_w_output(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, wpos, options);
fprintf("For inital guess weight: %0.5f: \n" + ...
    "Cost: %0.5f\n" + ...
    "Leg Max Height: %0.5f\n" + ...
    "Leg Min Height: %0.5f\n" + ...
    "Total positions: %i\n\n", wpos, out.cost, out.max_foot_height, out.min_foot_height, out.total_pos);
tableData = [tableData, Lw];
tableNames = [tableNames, "M5 We5"];
%% Mesh 9 Weight 0
wpos = 0;
Lw = [0.029488316154413	0.028561306370451	0.080738295904085	0.115399759429754	0.033875752011128	0.078396195002513	0.637983820867074]';
out = objective_w_output(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, wpos, options);
fprintf("For inital guess weight: %0.5f: \n" + ...
    "Cost: %0.5f\n" + ...
    "Leg Max Height: %0.5f\n" + ...
    "Leg Min Height: %0.5f\n" + ...
    "Total positions: %i\n\n", wpos, out.cost, out.max_foot_height, out.min_foot_height, out.total_pos);
tableData = [tableData, Lw];
tableNames = [tableNames, "M9 W0"];
runs{1} = out;
%% Mesh 9 Weight e3
wpos = 1e3;
Lw = [0.029958135975113	0.027654798769514	0.098426041132323	0.093915432374652	0.033679489783353	0.098753723868983	0.628958200876016]';
out = objective_w_output(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, wpos, options);
fprintf("For inital guess weight: %0.5f: \n" + ...
    "Cost: %0.5f\n" + ...
    "Leg Max Height: %0.5f\n" + ...
    "Leg Min Height: %0.5f\n" + ...
    "Total positions: %i\n\n", wpos, out.cost, out.max_foot_height, out.min_foot_height, out.total_pos);
tableData = [tableData, Lw];
tableNames = [tableNames, "M9 We3"];
runs{2} = out;
%% Mesh 9 Weight e4
wpos = 1e4;
Lw = [0.029653208536189	0.027681107202912	0.108047995757369	0.133199996441003	0.037799761503349	0.082400794588943	0.633392245411795]';
Lw2 = [0.029964840164286	0.027764051138324	0.101459898079935	0.131616339441701	0.037731970159884	0.080320029247764	0.718887343640291]';
out = objective_w_output(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, wpos, options);
out2 = objective_w_output(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw2, L0, wpos, options);
fprintf("For inital guess weight: %0.5f: \n" + ...
    "Cost: %0.5f\n" + ...
    "Leg Max Height: %0.5f\n" + ...
    "Leg Min Height: %0.5f\n" + ...
    "Total positions: %i\n\n", wpos, out.cost, out.max_foot_height, out.min_foot_height, out.total_pos);
tableData = [tableData, Lw];
tableNames = [tableNames, "M9 We4"];
runs{3} = out;
runs{5} = out2;
%% Mesh 9 Weight e5
wpos = 1e5;
Lw = [0.029907839572973	0.036275482997917	0.114174228736394	0.133199999691187	0.025268355002723	0.090974933341417	0.706699211077370]';
out = objective_w_output(Q0, Phi_system, Phi_q, nu, gamma, time, nb, omegas, Lw, L0, wpos, options);
fprintf("For inital guess weight: %0.5f: \n" + ...
    "Cost: %0.5f\n" + ...
    "Leg Max Height: %0.5f\n" + ...
    "Leg Min Height: %0.5f\n" + ...
    "Total positions: %i\n\n", wpos, out.cost, out.max_foot_height, out.min_foot_height, out.total_pos);
tableData = [tableData, Lw];
tableNames = [tableNames, "M9 We5"];
runs{4} = out;
%% Get table
fig1 = uifigure("Position",[500 500 800 400]);
T = array2table(tableData, 'VariableNames', tableNames, 'RowNames', names);
t = uitable(fig1, ...
    'Data', T, ...
    'Position', [20 20 760 360], ...
    'FontSize', 12, ...
    'FontName', 'Segoe UI', ...
    'RowStriping', 'on', ...
    'ColumnSortable', true);

% Plot data
figure;
bar(T{:,:})
grid
legend(tableNames)

figure;
bar([table2array(T(:,end-1)), Lw2])
grid
legend(["Pre correction", "Post correction"])

%%

MaxH = zeros(1,5);
MinH = zeros(1,5);
cost = zeros(1,5);
total_pos = zeros(1,5);

for i = 1:length(runs)
    if abs(runs{i}.max_foot_height) < 2e3
        MaxH(i) = runs{i}.max_foot_height;
    end
    MinH(i) = runs{i}.min_foot_height;
    if abs(runs{i}.cost) < 2e3
        cost(i) = runs{i}.cost;
    end
    
    total_pos(i) = runs{i}.total_pos;
end

legends = ["Weight: 0", "Weight: e3", "Weight: e4", "Weight: e5", "Corrected"];

fig2 = figure;
tld = tiledlayout(fig2, "flow");
nexttile;
bar(MaxH(1:5))
title("Max Height")
ylabel("Meters")
xticklabels(legends)

nexttile
bar(MinH(1:5))
title("Min Height")
ylabel("Meters")
xticklabels(legends)

nexttile;
bar(cost(1:5));
title("Cost")
ylabel("Cost")
xticklabels(legends)

nexttile;
bar(total_pos(1:5))
title("Total Positions")
ylabel("Configurations")
xticklabels(legends)