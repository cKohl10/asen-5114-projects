function nyquist_exp(exp_lg_mag, exp_lg_phase, model_lg, freq)
%NYQUIST_EXP Creates a Nyquist plot from the empirical data and overlays
%analytical Nyquist plot with stability margins.
% Inputs:
%   exp_lg_mag     - Empirical negative loop gain magnitude (linear scale)
%   exp_lg_phase   - Empirical negative loop gain phase (radians)
%   model_lg       - Analytical transfer function of the negative loop gain
%   freq           - empirical frequency vector


    % Convert empirical polar data to complex form
    emp_data = exp_lg_mag .* exp(1j * exp_lg_phase);
    
    % Get Nyquist data
    [re, im, wout] = nyquist(model_lg);
    re = squeeze(re);
    im = squeeze(im);
    [GM,PM] = margin(model_lg);
    phase_plot = PM-180;
    PM = [cosd(phase_plot),sind(phase_plot)];

    % Plot Nyquist curve
    figure;
    % Plot empirical Nyquist data
    exp_plot = plot(real(emp_data), imag(emp_data), 'r', 'LineWidth', 2, 'DisplayName', 'Empirical');
    hold on;
    plot(real(emp_data), -imag(emp_data), 'r', 'LineWidth', 2, 'DisplayName', 'Empirical');
    anal_plot = plot(re, im, 'b-', 'LineWidth', 2);        % Main loop
    plot(re, -im, 'b-', 'LineWidth', 2);      % Mirror for real systems
    
    % Critical point
    plot(-1, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    text(-1.5, 0, 'Crit Pt', 'FontSize', 12, 'Color', 'r');

    % Plot phase margin point (gain = 1 crossover)
    if exist('PM', 'var')
        plot(PM(1), PM(2), 'go', 'MarkerSize', 8, 'LineWidth', 2);
        text(PM(1), PM(2)-0.15, 'PM', 'Color', 'g', 'FontSize', 12);
    end

    % Plot gain margin point (phase = -180 crossover)
    if exist('GM', 'var')
        plot(-1/GM, 0, 'mo', 'MarkerSize', 8, 'LineWidth', 2);
        text(-1/GM, 0.1, '1/GM', 'Color', 'm', 'FontSize', 12);
    end

    theta = linspace(0, 2*pi, 500);
    x_unit = cos(theta);
    y_unit = sin(theta);
    
    plot(x_unit, y_unit, 'k--', 'LineWidth', 1.5);  
    axis equal;
    grid on;
    legend([exp_plot, anal_plot], 'Empirical', 'Analytical')
    xlabel('Real Axis', 'FontSize', 14);
    ylabel('Imaginary Axis', 'FontSize', 14);
    title('Loop Gain Nyquist Plot with Stability Margins', 'FontSize', 16);
    ax = gca;
    ax.FontSize = 12;
    xlim padded;
    ylim padded;
end

