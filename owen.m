% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig
clc;
clear;
close all;

% Read the data from the file
data = table2array(readtable('data/Spacecraft_spin_module_frequency_response_data.xlsx'));
freq = data(:,1);
amp = data(:,2);
phase = data(:,3);


% Estimate the transfer function
s = tf('s');
K = -1.5;
H_first = K/(s^2);
[mag_first,phase_first,w_first] = bode(H_first,freq);
mag_first = squeeze(mag_first);
phase_first = squeeze(phase_first);

% Plot the bode plot
figure;
hold on;
semilogx(freq, db(amp));
semilogx(w_first, db(mag_first));
title('Bode Plot of Spacecraft Spin Module Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
grid on;

% Plot the phase plot
figure;
semilogx(freq, rad2deg(phase));
title('Phase Plot of Spacecraft Spin Module Frequency Response');
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
grid on;


