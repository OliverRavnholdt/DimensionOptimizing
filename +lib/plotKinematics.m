% Takes in solved q, q dot and q double dot with the time and number of
% bodies
% Plots the position, velocity and acceleration for x, y and phi in a 3 x 3
% grid.

function fig1 = plotKinematics(q, qdot, qddot, t, nb)

% Setup plot
fig1 = figure;
tld = tiledlayout(fig1);
title(tld, "Kinematic analysis for one omega pair")

values = nb*3;
plotdata = {q, qdot, qddot};

pos_titles = ["X-Position", "Y-Position", "Angle"];
pos_ylabels = ["Position [m]", "Position [m]", "Angle [rad]"];

vel_titles = ["X-Velocity", "Y-Velocity", "Angular Velocity"];
vel_ylabels = ["Velocity [m/s]", "Velocity [m/s]", "Angular Velocity [rad/s]"];

acc_titles = ["X-Acceleration", "Y-Acceleration", "Angular Acceleration"];
acc_ylabels = ["Acceleration [m/s^2]", "Acceleration [m/s^2]", "Angular Acceleration [rad/s^2]"];

plot_titles = {pos_titles, vel_titles, acc_titles};
plot_ylabels = {pos_ylabels, vel_ylabels, acc_ylabels};

% Setup legends
legends = strings(1,nb);
for i = 1:nb
    legends(i) = sprintf("Body %i", i);
end

% Plot Pos, Vel, Acc
for i = 1:3
    % Plot X, Y, Angle
    for j = 1:3
        nexttile;
        title(plot_titles{i}(j));
        ylabel(plot_ylabels{i}(j));
        xlabel("Time [s]");
        
        % Plate data for each body
        for k = j:3:values
            hold on;

            if i == 1 && mod(k, 3) == 0
                plot(t, mod(plotdata{i}(k,:), 2*pi))
                continue;
            end
            plot(t, plotdata{i}(k,:))
        end
        legend(legends);
    end

end
end