function StepResponsePlot(sys, title_str)
% Calculate the step response
t = linspace(0, 10, 1000);
x0 = [0,0,0,0,0.1,0.1,0.1,0.1];
[y, t] = initial(sys, x0, t);
x_error = y(:,5:end);

% Figure setup
figure;
set(gcf, 'Position', [100, 100, 600, 800])
for i = 1:4
    subplot(4,1,i);
    plot(t, x_error(:,i), 'b-', 'LineWidth', 2); hold on;
    
    % Red line at 3 seconds
    xline(3, 'r--', 'LineWidth', 1.5, 'DisplayName', '3 sec marker');

    yline(0.05, 'k--', 'LineWidth', 1.5, 'DisplayName', '+5%')
    yline(-0.05, 'k--', 'LineWidth', 1.5, 'DisplayName', '-5%')

    xlabel('Time (s)', 'FontSize', 12);
    ylabel(['x_{error', num2str(i), '}(t)'], 'FontSize', 12);
    title(['State Error x_{', num2str(i), '}(t)'], 'FontSize', 14);
    ylim([-0.1,0.1])
    grid on;
    legend show;
end

sgtitle(['State Estimation Error (x\_error) - ' title_str], 'FontSize', 16);

end
