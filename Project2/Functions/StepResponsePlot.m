function StepResponsePlot(sys, title_str)
    [y, t] = step(sys);

    figure;
    set(gcf, 'Position', [100, 100, 500, 400]);

    plot(t, y, 'b-', 'LineWidth', 2); hold on;
    yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Steady-State (1)');
    xlabel('Time (s)', 'FontSize', 14);
    ylabel('Angle [rad]', 'FontSize', 14);
    title(['Step Response - ', title_str], 'FontSize', 16);
    grid on;

    % Annotate Overshoot
    % [peak_val, peak_idx] = max(y);
    % plot(t(peak_idx), peak_val, 'ro', 'MarkerSize', 8, 'DisplayName', 'Overshoot');
    % text(t(peak_idx), peak_val + 0.05, 'Overshoot', 'Color', 'r', 'FontSize', 12);

    % legend show;
end