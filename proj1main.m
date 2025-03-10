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

% Plant transfer function
DC_gain = 0.01;
p1 = 0.5; % rad/sec
p2 = 1; % rad/sec
s = tf('s');
P = (DC_gain)/(s^2*(s/p1+1)^2*(s/p2+1)^2);
[mag, phase, wout] = bode(P, freq);
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
semilogx(freq, rad2deg(phase_P0), '--');
title('Phase Plot');
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
xlim([freq(1), freq(end)]);
grid on;
legend('P', 'P_0');
saveas(gcf, 'figs/phase_plot.png');

