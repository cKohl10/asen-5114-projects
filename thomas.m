% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig4
clc;
clear;
close all;

% Read the data from the file
data = readmatrix('data\Spacecraft_spin_module_frequency_response_data.xlsx');
freq = data(:,1);
amp = data(:,2);
amp_db = db(amp) - 20*log10(freq);
phase = data(:,3) - pi/2;

% Create the analytical transfer function
s = tf('s');
K = db2mag(-1.5);
H_first = K / s^2;
[mag_1, phase_1, wout_1] = bode(H_first);


% Plot the bode plot
figure;
semilogx(freq, amp_db);
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


