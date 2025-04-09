%% Problem 2: Nyquist and Bode Plots Analysis
% This script plots the analytical and experiment nyquist and bode plots
% with various gain values for problem 2
% Author: Thomas Dunnington
close all; clear; clc;


%% Examine the analytical and empirical transfer functions with Nyquist and Bode
% Read the data from the file
data = readmatrix('data\Spacecraft_spin_module_frequency_response_data.xlsx');
freq_exp = data(:,1)*2*pi;
amp = data(:,2);
amp_db_P = 20*log10(amp) - 20*log10(freq_exp);
phase_P = data(:,3) - pi/2;

% Define the analytical transfer function
s = tf('s');
DC_gain = 10^(-15/20);
pole_1 = 0.35;
pole_2 = 0.35;
omega_ar1 = 4.601;   % anti-resonance frequency
omega_r1 = 8.347;    % resonance frequency
zeta_z = 0.015;
zeta_p = 0.035;
G_resonance = (s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2);
G_required = DC_gain * G_resonance * 1 / (s + pole_1) * 1 / (s + pole_2);




% Bode Plot of the natural system
[mag_model, phase_model, wout_model] = bode(G_required, {min(freq_exp), max(freq_exp)});
phase_model = squeeze(phase_model);
mag_model = squeeze(mag_model);

% Increase the gain
[mag_model_1, phase_model_1, wout_model_1] = bode(5*G_required, {min(freq_exp), max(freq_exp)});
phase_model_1 = squeeze(phase_model_1);
mag_model_1 = squeeze(mag_model_1);

[mag_model_2, phase_model_2, wout_model_2] = bode(10*G_required, {min(freq_exp), max(freq_exp)});
phase_model_2 = squeeze(phase_model_2);
mag_model_2 = squeeze(mag_model_2);

[mag_model_3, phase_model_3, wout_model_3] = bode(20*G_required, {min(freq_exp), max(freq_exp)});
phase_model_3 = squeeze(phase_model_3);
mag_model_3 = squeeze(mag_model_3);

[mag_model_4, phase_model_4, wout_model_4] = bode(100*G_required, {min(freq_exp), max(freq_exp)});
phase_model_4 = squeeze(phase_model_4);
mag_model_4 = squeeze(mag_model_4);


%% Bode
% Plot the experimental and analytical bode plots with different gain
% values
figure();
subplot(2,1,1)
semilogx(wout_model, db(mag_model), 'linewidth', 1.5)
hold on
semilogx(wout_model_1, db(mag_model_1), 'linewidth', 1.5)
semilogx(wout_model_2, db(mag_model_2), 'linewidth', 1.5)
semilogx(wout_model_3, db(mag_model_3), 'linewidth', 1.5)
semilogx(wout_model_4, db(mag_model_4), 'linewidth', 1.5)

grid on
title('Magnitude')
ylabel('Magnitude (dB)')
legend('Plant', 'K = 5', 'K = 10', 'K = 20', 'K = 100')

subplot(2,1,2)
semilogx(wout_model, phase_model, 'linewidth', 1.5);
hold on
yline(-180, 'color', 'r', 'linestyle', ':', 'LineWidth', 1.5)

grid on
title('Phase')
xlabel('Frequency (rad/s)')
ylabel('Phase (deg)')
sgtitle('Model Bode Plot')

%% Empirical Bode
figure();
subplot(2,1,1)
semilogx(freq_exp, amp_db_P, 'linewidth', 1.5)
hold on
semilogx(freq_exp, amp_db_P + db(5), 'linewidth', 1.5)
semilogx(freq_exp, amp_db_P + db(10), 'linewidth', 1.5)
semilogx(freq_exp, amp_db_P + db(20), 'linewidth', 1.5)
semilogx(freq_exp, amp_db_P + db(100), 'linewidth', 1.5)

grid on
title('Magnitude')
ylabel('Magnitude (dB)')
legend('Plant', 'K = 5', 'K = 10', 'K = 20', 'K = 100')

subplot(2,1,2)
semilogx(freq_exp, rad2deg(phase_P), 'linewidth', 1.5);
hold on
yline(-180, 'color', 'r', 'linestyle', ':', 'LineWidth', 1.5)

grid on
title('Phase')
xlabel('Frequency (rad/s)')
ylabel('Phase (deg)')
sgtitle('Empirical Bode Plot')


%% Large gain Bode plot
[mag_lg, phase_lg, wout_lg] = bode(5000*G_required, {min(freq_exp), max(freq_exp)});
mag_lg = squeeze(mag_lg);
phase_lg = squeeze(phase_lg);


% Plot the phase margin
gain_crossover_ind = find(diff(sign(db(mag_lg))) < 0 | diff(sign(db(mag_lg))) == 2);

% Plot the model bode plot
figure;
set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
subplot(2,1,1)
semilogx(wout_lg, db(mag_lg), 'b', 'linewidth', 2);
hold on;
yline(0, 'color', 'r', 'linestyle', ':', 'linewidth', 1.5)
xline(wout_lg(gain_crossover_ind), 'color', 'r', 'LineWidth', 1.5, 'linestyle', '--');
title('Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([wout_lg(1), wout_lg(end)]);
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(wout_lg, phase_lg, 'b', 'linewidth', 2);
hold on;
h1 = xline(wout_lg(gain_crossover_ind), 'color', 'r', 'LineWidth', 1.5, 'linestyle', '--');
yline(-180, 'color', 'r', 'linestyle', ':', 'linewidth', 1.5)
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([wout_lg(1), wout_lg(end)]);
legend([h1], 'Phase Margin Frequencies', 'location', 'nw')
sgtitle('Bode Plot')
grid on;


%% Nyquist
% Analytical and Empirical nyquist plots for the plant
[re,im,wout] = nyquist(G_required, {min(freq_exp), max(freq_exp)}); 
re = squeeze(re);
im = squeeze(im);
Gjw_exp = nyquist_exp(freq_exp, amp_db_P, phase_P, 0);

% Create Nyquist plot
figure();
hold on;
grid on;

% Experimental
p1 = plot(real(Gjw_exp), imag(Gjw_exp), 'b-', 'LineWidth', 1.5);
plot(real(Gjw_exp), -imag(Gjw_exp), 'b-', 'LineWidth', 1.5);

% Analytical
p2 = plot(re, im, 'r', 'LineWidth', 1.5);
plot(re, -im, 'r', 'LineWidth', 1.5);

% Formatting
axis equal;
xlabel('Real Axis', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Imaginary Axis', 'FontSize', 12, 'FontWeight', 'bold');
title('Experimental and Model Nyquist', 'FontSize', 14, 'FontWeight', 'bold');
legend([p1, p2], {'Experimental', 'Analytical'}, 'FontSize', 11);
set(gca, 'FontSize', 11)

% Grid Style
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);


%% Analytical Nyquist Plots with different gain values
figure();
nyquist(G_required, {min(freq_exp), max(freq_exp)}); 
hold on
% nyquist(5*G_required, {min(freq_exp), max(freq_exp)}); 
% nyquist(10*G_required, {min(freq_exp), max(freq_exp)});
% nyquist(20*G_required, {min(freq_exp), max(freq_exp)});
nyquist(5000*G_required, {min(freq_exp), max(freq_exp)}); 

axis equal
legend('Plant', 'K = 5000')
%legend('Plant', 'K = 5', 'K = 10', 'K = 20')


%% Experimental Nyquist Plots with Different Gain Values
% Calculate complex response for different magnitude values as a result of
% the proportional control
Gjw_exp_1 = nyquist_exp(freq_exp, amp_db_P, phase_P, 0);
Gjw_exp_2 = nyquist_exp(freq_exp, amp_db_P + db(5), phase_P, 0);
Gjw_exp_3 = nyquist_exp(freq_exp, amp_db_P + db(10), phase_P, 0);
Gjw_exp_4 = nyquist_exp(freq_exp, amp_db_P + db(20), phase_P, 0);

% Colors
default_colors = [
    0, 0.4470, 0.7410;  % Blue
    0.8500, 0.3250, 0.0980;  % Orange
    0.9290, 0.6940, 0.1250;  % Yellow
    0.4940, 0.1840, 0.5560;  % Purple
    0.4660, 0.6740, 0.1880;  % Green
    0.3010, 0.7450, 0.9330;  % Light Blue (Cyan)
    0.6350, 0.0780, 0.1840   % Dark Red
];

% Create Nyquist plot
figure();
hold on;
grid on;

% Experimental
p1 = plot(real(Gjw_exp_1), imag(Gjw_exp_1), 'color', default_colors(1,:), 'LineWidth', 1.5);
plot(real(Gjw_exp_1), -imag(Gjw_exp_1), 'color', default_colors(1,:), 'LineWidth', 1.5);

p2 = plot(real(Gjw_exp_2), imag(Gjw_exp_2), 'color', default_colors(2,:), 'LineWidth', 1.5);
plot(real(Gjw_exp_2), -imag(Gjw_exp_2), 'color', default_colors(2,:), 'LineWidth', 1.5);

p3 = plot(real(Gjw_exp_3), imag(Gjw_exp_3), 'color', default_colors(3,:), 'LineWidth', 1.5);
plot(real(Gjw_exp_3), -imag(Gjw_exp_3), 'color', default_colors(3,:), 'LineWidth', 1.5);

p4 = plot(real(Gjw_exp_4), imag(Gjw_exp_4), 'color', default_colors(4,:), 'LineWidth', 1.5);
plot(real(Gjw_exp_4), -imag(Gjw_exp_4), 'color', default_colors(4,:), 'LineWidth', 1.5);

scatter(-1, 0, 'Marker' , '+', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'SizeData', 100, 'LineWidth', 2);

% Formatting
axis equal;
xlabel('Real Axis', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Imaginary Axis', 'FontSize', 12, 'FontWeight', 'bold');
title('Experimental Nyquist', 'FontSize', 14, 'FontWeight', 'bold');
legend([p1, p2, p3, p4], {'Plant', 'K = 5', 'K = 10', 'K = 20'}, 'FontSize', 11);
set(gca, 'FontSize', 11)

% Grid Style
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);



% Just the plant
figure();
nyquist(G_required); 
axis equal

