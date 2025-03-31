% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig4
clc;
clear;
close all;

% Read the data from the file
data = readmatrix('data/Spacecraft_spin_module_frequency_response_data.xlsx');
freq = data(:,1)*2*pi;
amp = data(:,2);
amp_db_P = 20*log10(amp) - 20*log10(freq);
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
G = DC_gain * G_resonance * 1 / (s + pole_1) * 1 / (s + pole_2);

[mag, phase, wout] = bode(G, freq);
phase = squeeze(phase);
mag = squeeze(mag);

% Determine Open loop Poles and Zeros
OLPoles = pole(G);
OLZeros = zero(G);

% Plot the bode plot
% figure;
% semilogx(freq, amp_db_P);
% hold on;
% semilogx(wout, db(mag), 'r');
% title('Bode Plot');
% xlabel('Frequency (Hz)');
% ylabel('Amplitude (dB)');
% xlim([freq(1), freq(end)]);
% legend('Empirical', 'Analytical');
% grid on;
% % saveas(gcf, 'figs/bode_plot.png');
% 
% % Plot the phase plot
% figure;
% semilogx(freq, rad2deg(phase_P));
% hold on;
% semilogx(freq, phase, 'color', 'r');
% title('Phase Plot');
% xlabel('Frequency (rad/s)');
% ylabel('Phase (deg)');
% xlim([freq(1), freq(end)]);
% grid on;
% legend('Experimental', 'Analytical');
% % saveas(gcf, 'figs/phase_plot.png');


% %%%%% Problem 2 %%%%%
% Simple controller - plot all K values on one Nyquist plot
% K = [-1000, -100, -10, -1, -0.1 ,0.1, 1, 10, 100, 1000];
% figure;
% hold on;
% legendEntries = cell(1, length(K));
% 
% for i = 1:length(K)
%     C = K(i);
%     nyquist(C*G);
%     legendEntries{i} = sprintf('K = %g', K(i));
% end
% 
% title('Nyquist Plot - Multiple Gain Values');
% xlabel('Real Axis');
% ylabel('Imaginary Axis');
% legend(legendEntries);
% grid on;
% 
% 
% % Simple controller - plot all K values on one Bode plot
% figure;
% hold on;
% legendEntries = cell(1, length(K));
% 
% for i = 1:length(K)  
%     C = K(i);
%     clSys = feedback(C*G, 1);  % closed-loop transfer function
%     polesCL(:,i) = pole(clSys);
%     bode(C*G);
%     legendEntries{i} = sprintf('K = %g', K(i));
% end
% 
% title('Bode Plot - Multiple Gain Values');
% legend(legendEntries);
% grid on;
% % saveas(gcf, 'figs/nyquist_plot_combined.png');
% 

%%%%% Problem 3 %%%%%
% Design a compensator to meet the following design requirements:
% Phase Margin of > 40 deg
% Gain margin > 10 dB
% CL Bandwidth ~ 1 Hz

% Define Compensator Properties
z1 = 0.01;   
p1 = 6*pi; 
z2 = 1;
p2 = 6*pi;

% Compensator Gain (Tuned for Stability)
K = 8000;
zeta_z = 0.01;
zeta_p = 0.06;

% Combined Notch Filter to Damp Resonance and Anti-Resonance
notch = (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2)/(s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2);
C = K * (s + z1)/(s + p1) * (s + z2)/(s + p2) * 1/s * notch;
% Calculate Gain and Phase Margin 
Lg = -C*G;
[GM,PM] = margin(-Lg); 
GM = 20*log10(GM);
Margin = [GM,PM]';

% Determine Closed Loop Poles and Zeros
clSys = feedback(-Lg, 1);  
CLPoles = pole(clSys);
CLZeros = zero(clSys);
bw = bandwidth(clSys);

% Evaluate Tracking 
w = logspace(-1, 2.5, 1000);  
for i  = 1:length(w)
    Lg_eval(i) = evalfr(Lg, w(i));
    Ts_eval(i) = 1/abs(1-Lg_eval(i)); 
end

% Display the margins and bandwidth 
disp(['Gain Margin: ', num2str(GM), 'dB']);
disp(['Phase Margin: ', num2str(PM), 'deg']);
disp(['Closed-loop bandwidth: ', num2str(bw), ' rad/s']);

% Create a Bode Plot for Loop Gain
figure;
hold on;
bode(-Lg);
title('Loop Gain Bode Plot with Compensator');
grid on;
% saveas(gcf, 'figs/P3_BodePlot.png');

% Create a Bode Plot for Loop Gain
figure;
hold on;
bode(clSys);
title('Closed Loop Bode Plot with Compensator');
grid on;
% saveas(gcf, 'figs/P3_BodePlot.png');

% Create a Nyquist Plot
figure;
hold on;
nyquist(-Lg);
title('Loop Gain Nyquist Plot')
grid on;

% Plot Poles and Zeros 
figure;
hold on;
plot(real(CLZeros), imag(CLZeros), 'ro', 'MarkerSize', 8);
plot(real(CLPoles), imag(CLPoles), 'rx', 'MarkerSize', 10);
plot(real(OLZeros), imag(OLZeros), 'bo', 'MarkerSize', 8);
plot(real(OLPoles), imag(OLPoles), 'bx', 'MarkerSize', 10);

xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Pole-Zero Plot');
grid on;
axis equal; 
legend('Closed Loop Zeros', 'Closed Loop Poles','Open Loop Zero','Open Loop Poles');
hold off; 


% Plot Tracking Error
figure;
semilogx(w, Ts_eval*100, 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('Tracking Error (%)');
title('Closed-Loop Tracking Error vs Frequency');
grid on;
% yline(5, '--r', '5% Error Threshold');
ylim([0 100]);
