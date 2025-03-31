%% Compensation Design
% This script designs a compensator to achieve design objectives
close all; clear; clc;

%% Plant
% Define the analytical transfer function
s = tf('s');
DC_gain = 10^(-15/20);
pole_1 = 0.8;
pole_2 = 0.01;
omega_ar1 = 4.601;  
omega_r1 = 8.347;  
zeta_z = 0.015;
zeta_p = 0.03;
G_resonance = (s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2);
G = DC_gain * G_resonance * 1 / (s + pole_1) * 1 / (s + pole_2);

% Empirical Frequency Ranges
data = readmatrix('data\Spacecraft_spin_module_frequency_response_data.xlsx');
freq_exp = data(:,1)*2*pi;
w_min = min(freq_exp);
w_max = max(freq_exp);

%% Compensation
%%%% Notch Filter
wn = omega_ar1;        % Notch frequency
zeta = 0.03;     
zeta_prime = 0.2; 

% Notch Filter Transfer Function
notch_filter = (s^2 + 2*zeta_prime*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta*omega_ar1*s + omega_ar1^2);

% Visualize the notch
% figure();
% bode(notch_filter);


% Note: Try to achieve control objectives without using a notch filter to
% damp out the anti-resonance first (i.e. just use lead and constant gain)


%%%% Lead compensator
K = 1000;
a1 = 2;
b1= 5.5;
C1 = K * (s + a1) / (s + b1);

% Second compensator
a2 = 7;
b2 = 25;
C2 = b2/a2 * (s + a2) / (s+ b2);

% Controller
C = C1*C2;

% Negative Loop Gain
Lg = C*G;

% Compute the frequency response
[mag_lg, phase_lg, wout_lg] = bode(Lg, {w_min, w_max});
phase_lg = squeeze(phase_lg);
mag_lg = squeeze(mag_lg);

% Calculate the closed loop response
T = feedback(Lg, 1);
[mag_cl, phase_cl, wout_cl] = bode(T, {w_min, w_max});
phase_cl = squeeze(phase_cl);
mag_cl = squeeze(mag_cl);

% Calculate the closed loop bandwidth
closed_loop_bandwidth = bandwidth(T);

% Calculate the open loop stability margins
% Phase margin
gain_crossover_ind = find(diff(sign(db(mag_lg))) < 0 | diff(sign(db(mag_lg))) == 2);

% Gain Margin
phase_crossover_ind = find(diff(sign(phase_lg + 180)) < 0 | diff(sign(phase_lg + 180)) == 2);

%% Loop Gain Bode Plot
% Plot the model bode plot
figure;
set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
subplot(2,1,1)
semilogx(wout_lg, db(mag_lg), 'b', 'linewidth', 2);
hold on;
if(~isempty(phase_crossover_ind))
    for i = 1:length(phase_crossover_ind)
        if abs(db(mag_lg(phase_crossover_ind(i)))) > 10
            line([wout_lg(phase_crossover_ind(i)), wout_lg(phase_crossover_ind(i))], [0, db(mag_lg(phase_crossover_ind(i)))], 'Color', 'g', 'LineWidth', 1.5);
        else
            line([wout_lg(phase_crossover_ind(i)), wout_lg(phase_crossover_ind(i))], [0, db(mag_lg(phase_crossover_ind(i)))], 'Color', 'r', 'LineWidth', 1.5);
        end
    end
end
yline(0, 'color', 'r', 'linestyle', ':', 'linewidth', 1.5)
title('Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([wout_lg(1), wout_lg(end)]);
legend('Compensated Response');
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(wout_lg, phase_lg, 'linewidth', 2);
hold on;
for i = 1:length(gain_crossover_ind)
    % Plot green if good phase margin, otherwise red
    if phase_lg(gain_crossover_ind(i)) >= -140
        line([wout_lg(gain_crossover_ind(i)), wout_lg(gain_crossover_ind(i))], [-180, phase_lg(gain_crossover_ind(i))], 'Color', 'g', 'LineWidth', 1.5);
    else
        line([wout_lg(gain_crossover_ind(i)), wout_lg(gain_crossover_ind(i))], [-180, phase_lg(gain_crossover_ind(i))], 'Color', 'r', 'LineWidth', 1.5);
    end
end
yline(-180, 'color', 'r', 'linestyle', ':', 'linewidth', 1.5)
yline(-140, 'color', 'g', 'linestyle', ':', 'linewidth', 1.5)
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([wout_lg(1), wout_lg(end)]);
sgtitle('Loop Gain Bode Plot')
grid on;


%% Closed Loop Bode Plot
figure;
set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
subplot(2,1,1)
semilogx(wout_cl, db(mag_cl), 'color', 'g', 'linewidth', 2);
hold on;
xline(closed_loop_bandwidth, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--', 'label', closed_loop_bandwidth)
yline(-3, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
title('Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([wout_cl(1), wout_cl(end)]);
legend('Compensated Response', 'Bandwidth');
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(wout_cl, phase_cl, 'color', 'g', 'linewidth', 2);
hold on;
% xline(2*pi, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([wout_cl(1), wout_cl(end)]);
sgtitle('Closed Loop Bode Plot')
grid on;


%% Nyquist Plot
% Create a nyquist plot of the model
figure();
nyquist(Lg, {w_min, w_max});
axis equal;

