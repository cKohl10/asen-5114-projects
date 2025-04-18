function control_plot(out, fig_size, fig_save, title_str)
    f = figure();
    f.Position = fig_size;
    subplot(2, 1, 1);
    plot(out.input.Time, out.input.Data, 'r', 'LineWidth', 2, 'DisplayName', 'Reference Angle (deg)');
    hold on;
    plot(out.theta.Time, out.theta.Data, 'b', 'LineWidth', 2, 'DisplayName', 'Response (deg)');
    title(sprintf('Response - %s', title_str));
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    grid on;
    legend('show');

    subplot(2, 1, 2);
    plot(out.u.Time, out.u.Data, 'r', 'LineWidth', 2, 'DisplayName', 'Torque (mNm)');
    title(sprintf('Controls - %s', title_str));
    xlabel('Time (s)');
    ylabel('Torque (mNm)');
    grid on;
    legend('show');
    saveas(f, fig_save + "/" + title_str + ".png");
end
