% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig4
clc;
clear;
close all;

% Read the data from the file
data = table2array(readtable('data\Spacecraft_spin_module_frequency_response_data.xlsx'));
freq = data(:,1);
amp = data(:,2);
amp_db_P = 20*log10(amp) - 20*log10(freq);
amp_db_P0 = 20*log10(amp);
phase_P = data(:,3) - pi/2;
phase_P0 = data(:,3);

% Transfer function with anti resonance, resonance pair
s = tf('s');
DC_gain = 0.041;
omega_ar1 = 0.73;   % anti-resonance frequency
omega_r1 = 1.32;    % resonance frequency
omega_ar2 = 4.9;
omega_r2 = 6.1;
zeta_z = 0.015;
zeta_p = 0.03;
zeta_z2 = 0.03;
zeta_p2 = 0.034;
G = DC_gain * (1/s^2) * (s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2) * (s^2 + 2*zeta_z2*omega_ar2*s + omega_ar2^2) / (s^2 + 2*zeta_p2*omega_r2*s + omega_r2^2);
[mag, phase, wout] = bode(G, freq);
phase = squeeze(phase);
mag = squeeze(mag);


% Plot the bode plot
figure;
semilogx(freq, amp_db_P);
hold on;
semilogx(freq, amp_db_P0, '--');
semilogx(wout, db(mag), 'r');
title('Bode Plot');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
xlim([freq(1), freq(end)]);
legend('P', 'P_0', 'P(s)');
grid on;
saveas(gcf, 'figs/bode_plot.png');

% Plot the phase plot
figure;
semilogx(freq, rad2deg(phase_P));
hold on;
semilogx(freq, phase, '--');
title('Phase Plot');
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
xlim([freq(1), freq(end)]);
grid on;
legend('Experimental', 'Analytical');
saveas(gcf, 'figs/phase_plot.png');


%% Nyquist plotting
nyquist_exp(freq, amp_db_P, phase);