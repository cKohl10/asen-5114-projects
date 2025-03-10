% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig4
clc;
clear;
close all;

% Read the data from the file
data = table2array(readtable('data\Spacecraft_spin_module_frequency_response_data.xlsx'));
freq = data(:,1);
amp = data(:,2);
amp_db = 20*log10(amp) - 20*log10(freq);
phase = data(:,3) - pi/2;

% Plot the bode plot
figure;
semilogx(freq, amp_db);
title('Bode Plot of Spacecraft Spin Module Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
grid on;
saveas(gcf, 'figs/bode_plot.png');

% Plot the phase plot
figure;
semilogx(freq, rad2deg(phase));
title('Phase Plot of Spacecraft Spin Module Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on;
saveas(gcf, 'figs/phase_plot.png');

