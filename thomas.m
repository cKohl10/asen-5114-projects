% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig4
clc;
clear;
close all;

% Read the data from the file
data = readmatrix('data\Spacecraft_spin_module_frequency_response_data.xlsx');
freq_exp = data(:,1)*2*pi;
amp = data(:,2);
amp_db_P = 20*log10(amp) - 20*log10(freq_exp);
amp_db_P0 = 20*log10(amp);
phase_P = data(:,3) - pi/2;
phase_P0 = data(:,3);

% Transfer function with anti resonance, resonance pair
s = tf('s');
DC_gain = 10^(-15/20);
pole_1 = 0.8;
pole_2 = 0.01;
omega_ar1 = 4.601;   % anti-resonance frequency
omega_r1 = 8.347;    % resonance frequency
zeta_z = 0.015;
zeta_p = 0.03;
G_resonance = (s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2);
G_required = DC_gain * G_resonance * 1 / (s + pole_1) * 1 / (s + pole_2);


% Intermediate model
DC_gain = 10^(-15/20);
pole_1 = 0.3;
pole_2 = 0.4;
omega_ar1 = 4.601;   % anti-resonance frequency
omega_r1 = 8.347;    % resonance frequency
zeta_z = 0.015;
zeta_p = 0.035;
G_resonance = (s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2);
G_required = DC_gain * G_resonance * 1 / (s + pole_1) * 1 / (s + pole_2);



% Bode
figure();
bode(G_required)

% Compute the frequency response
[mag, phase, wout] = bode(G_required, {min(freq_exp), max(freq_exp)});
phase = squeeze(phase);
mag = squeeze(mag);


%% Only Empirical Bode
% Plot the bode plot
figure;
set(gcf, 'Position', [100, 100, 600, 500]); % Resize figure window
subplot(2,1,1)
semilogx(freq_exp, amp_db_P, 'linewidth', 2);
title('Magnitude');
ylabel('Amplitude (dB)');
xlim([freq_exp(1), freq_exp(end)]);
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(freq_exp, rad2deg(phase_P), 'linewidth', 2);
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([freq_exp(1), freq_exp(end)]);
sgtitle('Bode Plot')
grid on;


%% Full Bode Plot

% Plot the bode plot
figure;
set(gcf, 'Position', [100, 100, 600, 500]); % Resize figure window
subplot(2,1,1)
semilogx(freq_exp, amp_db_P, 'linewidth', 2);
hold on;
semilogx(wout, db(mag), 'r', 'linewidth', 2);
title('Magnitude');
ylabel('Amplitude (dB)');
xlim([freq_exp(1), freq_exp(end)]);
legend('Empirical', 'Analytical');
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(freq_exp, rad2deg(phase_P), 'linewidth', 2);
hold on;
semilogx(wout, phase, 'color', 'r', 'linewidth', 2);
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([freq_exp(1), freq_exp(end)]);
sgtitle('Bode Plot')
grid on;

%saveas(gcf, 'figs/bode_plot_full.png');


%% Nyquist plotting
%nyquist_exp(freq_exp, amp_db_P, phase_P, 1);


