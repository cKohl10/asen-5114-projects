function Gjw = nyquist_exp(freq, mag, phase, plotBool)
%NYQUIST_EXP Creates a Nyquist plot from the empirical data
% freq = rad/s
% mag = db(rad/s/V)
% phase = rad
% plotBool = true -> Creates the figure, otherwise just returns the
% complex response

% Convert magnitude to linear scale
mag_linear = 10.^(mag / 20);

% Convert phase to rad
phase_rad = phase;

% Compute complex response
Gjw = mag_linear .* exp(1j * phase_rad);

if (plotBool)
    % Nyquist plot
    figure;
    plot(real(Gjw), imag(Gjw), 'b', 'LineWidth', 1.5);
    hold on;
    plot(real(Gjw), -imag(Gjw), 'r--', 'LineWidth', 1.2); % Mirror for coming back on the Nyquist path
    grid on;
    xlabel('Re(G(j\omega))');
    ylabel('Im(G(j\omega))');
    title('Empirical Nyquist Plot');
    axis equal;
    
    % Add arrows for the direction
    quiver(real(Gjw(1:end-1)), imag(Gjw(1:end-1)), ...
       diff(real(Gjw)), diff(imag(Gjw)), 0, 'k');
    
    % Highlight the critical point
    % scatter(-1,0,'filled','marker','*', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'SizeData', 40)
    
    
    % plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
    % [~, idx] = min(abs(Gjw + 1)); % Closest point to -1 (stability check)
    % plot(real(Gjw(idx)), imag(Gjw(idx)), 'ro', 'MarkerFaceColor', 'r');
end
end

