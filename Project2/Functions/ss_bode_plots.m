function ss_bode_plots(A, B, C, D, K, F)
%SS_BODE_PLOTSS Takes as input a state space model, state feedback gain
%matrix K, and feed forward gain matrix F and plots the loop gain and
%closed loop bode plots
% Inputs:
%   A, B, C, D -> State space model
%   K -> Feedback gain matrix
%   F -> Feedforward gain matrix (Make F a scalar for a single input
%   system)
% Outputs:
%   A figure of the loog gain bode plot with stability margins
%   A figure of the closed loop bode plot with bandwidth

% Loop gain tf
s = tf('s');
Lg = tf(ss(A, B, K, 0));

% Closed loop tf
Acl = A - B*K;
sys_cl = ss(Acl, B, C, 0);
Cl = tf(sys_cl) * F;

% Torque input closed loop tf
In = inv(eye(length(K(:,1))) + Lg)*F;

% Calculate the frequency responses
[mag_lg, phase_lg, wout_lg] = bode(Lg);
phase_lg = squeeze(phase_lg);
mag_lg = squeeze(mag_lg);

[mag_cl, phase_cl, wout_cl] = bode(Cl);
phase_cl = squeeze(phase_cl);
mag_cl = squeeze(mag_cl);

[mag_in, phase_in, wout_in] = bode(In);
phase_in = squeeze(phase_in);
mag_in = squeeze(mag_in);

% Calculate the closed loop bandwidth
closed_loop_bandwidth = bandwidth(Cl,-3.005);

% Calculate the open loop stability margins
% Phase margin
gain_crossover_ind = find(diff(sign(db(mag_lg))) < 0 | diff(sign(db(mag_lg))) == 2);

% Gain Margin
phase_crossover_ind = find(diff(sign(phase_lg + 180)) < 0 | diff(sign(phase_lg + 180)) == 2 | diff(sign(phase_lg - 180)) < 0 | diff(sign(phase_lg - 180)) == 2);

% Report crossover frequencies and minimum margins
if ~isempty(gain_crossover_ind)
    min_pm = 360;
    min_idx = gain_crossover_ind(1);
    for i = 1:length(gain_crossover_ind)
        idx = gain_crossover_ind(i);
        current_pm = min([abs(phase_lg(idx) - 180), abs(phase_lg(idx) + 180)]);
        if (current_pm < min_pm)
            min_pm = current_pm;
            min_idx = idx;
        end
    end
    fprintf('Phase Margin: %.2f deg (at %.3f rad/s)\n', min_pm, wout_lg(min_idx)); 
else
    fprintf('No Gain Crossover Frequency found.\n');
end

if ~isempty(phase_crossover_ind)
    idx = phase_crossover_ind(1); % Report based on the first crossing
    current_gm = -db(mag_lg(idx)); % GM in dB
    fprintf('Gain Margin: %.2f dB (at %.3f rad/s)\n', current_gm, wout_lg(idx)); 
else
    fprintf('No Phase Crossover Frequency found (Gain Margin is potentially Infinite).\n');
end
fprintf('Closed Loop Bandwidth: %.3f rad/s\n', closed_loop_bandwidth);

%% Loop Gain Bode Plot
% Plot the model bode plot
figure;
set(gcf, 'Position', [100, 100, 500, 400]); % Resize figure window
subplot(2,1,1)
semilogx(wout_lg, db(mag_lg), 'b', 'linewidth', 2);
hold on;
yline(-10, 'color', 'g', 'linestyle', ':', 'linewidth', 1.5)
if(~isempty(phase_crossover_ind))
    for i = 1:length(phase_crossover_ind)
        if db(mag_lg(phase_crossover_ind(i))) < -10
            good_gain = line([wout_lg(phase_crossover_ind(i)), wout_lg(phase_crossover_ind(i))], [0, db(mag_lg(phase_crossover_ind(i)))], 'Color', 'g', 'LineWidth', 1.5);
        else
            bad_gain = line([wout_lg(phase_crossover_ind(i)), wout_lg(phase_crossover_ind(i))], [0, db(mag_lg(phase_crossover_ind(i)))], 'Color', 'r', 'LineWidth', 1.5);
        end
    end
end
yline(0, 'color', 'r', 'linestyle', ':', 'linewidth', 1.5)

% Add scatter points for crossovers
if ~isempty(gain_crossover_ind)
    scatter(wout_lg(gain_crossover_ind), zeros(size(gain_crossover_ind)), 60, 'm', 'filled', 'DisplayName', 'Gain Crossover');
end
if ~isempty(phase_crossover_ind)
    scatter(wout_lg(phase_crossover_ind), db(mag_lg(phase_crossover_ind)), 60, 'c', 'filled', 'DisplayName', 'Phase Crossover');
end

title('Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([wout_lg(1), wout_lg(end)]);
grid on;

% Add legend
try
    legend([good_gain, bad_gain], 'Good Gain Margin', 'Bad Gain Margin', 'location' ,'best');
catch
    try 
        legend([good_gain], 'Good Gain Margin', 'location' ,'best');
    catch
        try
            legend([bad_gain], 'Bad Gain Margin', 'location' ,'best');
        catch
            x = 1;
        end
    end
end

% Plot the phase plot
subplot(2,1,2)
semilogx(wout_lg, phase_lg, 'linewidth', 2);
hold on;
for i = 1:length(gain_crossover_ind)
    % Check if the phase margin is good or bad
    if (abs(phase_lg(gain_crossover_ind(i)) - 180) >= 40 && abs(phase_lg(gain_crossover_ind(i)) + 180) >= 40)
        if (abs(phase_lg(gain_crossover_ind(i)) - 180) < abs(phase_lg(gain_crossover_ind(i)) + 180))
            good_phase = line([wout_lg(gain_crossover_ind(i)), wout_lg(gain_crossover_ind(i))], [180, phase_lg(gain_crossover_ind(i))], 'Color', 'g', 'LineWidth', 1.5);
        else
            good_phase = line([wout_lg(gain_crossover_ind(i)), wout_lg(gain_crossover_ind(i))], [-180, phase_lg(gain_crossover_ind(i))], 'Color', 'g', 'LineWidth', 1.5);
        end
    else
        % Closest to 180 degrees
        if (abs(phase_lg(gain_crossover_ind(i)) - 180) < 40)
            bad_phase = line([wout_lg(gain_crossover_ind(i)), wout_lg(gain_crossover_ind(i))], [180, phase_lg(gain_crossover_ind(i))], 'Color', 'r', 'LineWidth', 1.5);
        else
            bad_phase = line([wout_lg(gain_crossover_ind(i)), wout_lg(gain_crossover_ind(i))], [-180, phase_lg(gain_crossover_ind(i))], 'Color', 'r', 'LineWidth', 1.5);
        end
    end
end
yline(-180, 'color', 'r', 'linestyle', ':', 'linewidth', 1.5)
yline(180, 'color', 'r', 'linestyle', ':', 'linewidth', 1.5)
yline(-140, 'color', 'g', 'linestyle', ':', 'linewidth', 1.5)

% Add scatter points for crossovers
if ~isempty(gain_crossover_ind)
    scatter(wout_lg(gain_crossover_ind), phase_lg(gain_crossover_ind), 60, 'm', 'filled', 'DisplayName', 'Gain Crossover');
end
if ~isempty(phase_crossover_ind)
    scatter(wout_lg(phase_crossover_ind), phase_lg(phase_crossover_ind), 60, 'c', 'filled', 'DisplayName', 'Phase Crossover');
end

title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([wout_lg(1), wout_lg(end)]);
sgtitle('Loop Gain Bode Plot')
grid on;

% Add legend
try
    legend([good_phase, bad_phase], 'Good Phase Margin', 'Bad Phase Margin', 'location' ,'best');
catch
    try 
        legend([good_phase], 'Good Phase Margin', 'location' ,'best');
    catch
        try
            legend([bad_phase], 'Bad Phase Margin', 'location' ,'best');
        catch
            x = 1;
        end
    end
end


%% Closed Loop Bode Plot
figure;
set(gcf, 'Position', [100, 100, 500, 400]); % Resize figure window
subplot(2,1,1)
semilogx(wout_cl, db(mag_cl), 'color', 'g', 'linewidth', 2);
hold on;
xline(2*pi, 'k--', 'LineWidth', 1, 'DisplayName', '1 Hz', 'LabelVerticalAlignment', 'bottom'); % Add line at 1 Hz with label
bandwidth_plot = xline(closed_loop_bandwidth, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--', 'label', closed_loop_bandwidth, 'LabelVerticalAlignment', 'bottom');
yline(-3, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
title('Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([wout_cl(1), wout_cl(end)]);
legend([bandwidth_plot], 'Bandwidth');
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(wout_cl, phase_cl, 'color', 'g', 'linewidth', 2);
hold on;
xline(2*pi, 'k--', 'LineWidth', 1, 'DisplayName', '1 Hz', 'Label', '1 Hz', 'LabelVerticalAlignment', 'bottom'); % Add line at 1 Hz with label
% xline(2*pi, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([wout_cl(1), wout_cl(end)]);
sgtitle('Closed Loop Bode Plot')
grid on;

%% Control Input Bode Plot
figure;
set(gcf, 'Position', [100, 100, 500, 400]); % Resize figure window
subplot(2,1,1)
semilogx(wout_in, db(mag_in), 'color', 'r', 'linewidth', 2);
hold on;
yline(20*log10(67), 'k--', 'LineWidth', 1, 'label', 'Saturation Limit', 'LabelHorizontalAlignment', 'left'); % Add line at 1 Hz with label
title('Magnitude (mN-m/rad)');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([wout_in(1), wout_in(end)]);
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(wout_in, phase_in, 'color', 'r', 'linewidth', 2);
hold on;
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([wout_in(1), wout_in(end)]);
sgtitle('Reference Input to Plant Input')
grid on;

end

