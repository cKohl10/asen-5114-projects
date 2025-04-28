function control_plot(out, fig_size, fig_save, title_str)
    f = figure();
    fs = 12;
    f.Position = fig_size;
    subplot(2, 1, 1);
    plot(out.input.Time, out.input.Data, 'r', 'LineWidth', 2, 'DisplayName', 'Reference Angle (deg)');
    set(gca, 'FontSize', fs)
    hold on;
    plot(out.theta.Time, out.theta.Data, 'b', 'LineWidth', 2, 'DisplayName', 'Response (deg)');
    set(gca, 'FontSize', fs)
    title(sprintf('Response - %s', title_str), 'FontSize', fs+2);
    xlabel('Time (s)', 'FontSize', fs);
    ylabel('Angle (deg)', 'FontSize', fs);
    grid on;
    legend('show');

    subplot(2, 1, 2);
    plot(out.u.Time, out.u.Data, 'r', 'LineWidth', 2, 'DisplayName', 'Torque (mNm)');
    set(gca, 'FontSize', fs)
    title(sprintf('Controls - %s', title_str), 'FontSize', fs+2);
    xlabel('Time (s)', 'FontSize', fs);
    ylabel('Torque (mNm)', 'FontSize', fs);
    grid on;
    legend('show');
    saveas(f, fig_save + "/" + title_str + ".png");
end
