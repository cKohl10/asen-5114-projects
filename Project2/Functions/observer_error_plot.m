function observer_error_plot(t, x, xhat, sim_title)
%OBSERVER_ERROR_PLOT Plots the true state, the observed state, and the
%state error for a given simulation.
% Inputs:
%   t -> time
%   x -> true state (p x 4)
%   xhat -> observed state (p x 4)
%   sim_title -> title string for the type of simulation

% Calculate the state error
xe = x - xhat;

% Plot formatting parameters
lw = 2;
fs = 12;

% State plot
figure('Name', 'Observed and True States', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
sgtitle(['Observed and True States: ', sim_title], 'FontSize', fs+2)

for i = 1:4
    subplot(4,1,i)
    plot(t, x(:,i), 'Color', 'b', 'LineWidth', lw)
    hold on
    plot(t, xhat(:,i), 'Color', 'r', 'LineWidth', lw, 'LineStyle', '--')
    ylabel(['x_', num2str(i)], 'FontSize', fs)
    grid on
    if i == 1
        legend('True', 'Estimated', 'Location', 'best')
    end
end
xlabel('Time (s)', 'FontSize', fs)

% State error plot
figure('Name', 'Observer Error', 'NumberTitle', 'off', 'Position', [100, 100, 800, 600]);
sgtitle(['Observer Error: ', sim_title], 'FontSize', fs+2)

for i = 1:4
    subplot(4,1,i)
    plot(t, xe(:,i), 'Color', 'r', 'LineWidth', lw)
    ylabel(['e_', num2str(i)], 'FontSize', fs)
    grid on
end
xlabel('Time (s)', 'FontSize', fs)

end
