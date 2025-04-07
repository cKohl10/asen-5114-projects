% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig
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
pole_1 = 0.3;
pole_2 = 0.4;
omega_ar1 = 4.601;   % anti-resonance frequency
omega_r1 = 8.347;    % resonance frequency
zeta_z = 0.015;
zeta_p = 0.035;
G_resonance = (s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2);
G = DC_gain * G_resonance * 1 / (s + pole_1) * 1 / (s + pole_2);


[mag, phase, wout] = bode(G, freq);
phase = squeeze(phase);
mag = squeeze(mag);

% Determine Open loop Poles and Zeros
OLPoles = pole(G);
OLZeros = zero(G);

% Bode Plot of the Open Loop
figure;
subplot(2,1,1)
semilogx(freq, amp_db_P, 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('|-Lg(j\omega)| (dB)');
title('Magnitutde of the Empirical Open-Loop');
grid on;
subplot(2,1,2)
semilogx(freq, phase_P, 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
title('Phase of Empirical Open-Loop Response');
grid on;

%%%%% Problem 3 %%%%%
% Design a compensator to meet the following design requirements:
% Phase Margin of > 40 deg
% Gain margin > 10 dB
% CL Bandwidth ~ 1 Hz

% Define Compensator
% Pole Cancellation: Cancel out the zero and pole with "approximate" cancellation
zeta_r = 2*0.035;
zeta_ar = 0.5*0.015;
pole_cancel = (s^2 + 2*zeta_r*omega_r1*s + omega_r1^2) / (s^2 + 2*zeta_ar*omega_ar1*s + omega_ar1^2);

% Proportional Gain
K = 25;

% Lead compensator
a1 = 1;
b1 = 10;
C1 = b1 / a1 * (s + a1) / (s + b1);

% Total Compensator
C = K*C1*pole_cancel;

% Negative Loop Gain
Lg = C*G;

% Calculate Gain and Phase Margin 
[GM,PM] = margin(Lg); 
GM = 20*log10(GM);
Margin = [GM,PM]';

% Determine Closed Loop Poles and Zeros
clSys = feedback(Lg, 1);  
CLPoles = pole(clSys);
CLZeros = zero(clSys);
bw = bandwidth(clSys);

% Evaluate Tracking 
w = logspace(-1, 2.5, 1000);  
for i  = 1:length(w)
    Lg_eval(i) = evalfr(-Lg, w(i));
    Ts_eval(i) = 1/abs(1-Lg_eval(i)); 
end

% Display the margins and bandwidth 
disp(['Gain Margin: ', num2str(GM), 'dB']);
disp(['Phase Margin: ', num2str(PM), 'deg']);
disp(['Closed-loop bandwidth: ', num2str(bw), ' rad/s']);

% Create a Bode Plot for Loop Gain
% figure;
% hold on;
% bode(Lg);
% title('Loop Gain Bode Plot with Compensator');
% grid on;
% % saveas(gcf, 'figs/P3_BodePlot.png');
% 
% % % Create Bode plot of compesators and LG
% % w = logspace(-1, 2, 1000);  % 0.1 to 100 rad/s
% % figure;
% % bode(C, w); 
% % hold on;
% % bode(Lg,w)
% % bode(G,w)
% % legend('C','LG','Plant');
% % title('Bode Plots of Compensator, LG, Plant');
% % grid on;
% 
% 
% % Create a Bode Plot for Loop Gain
% figure;
% hold on;
% bode(clSys);
% title('Closed Loop Bode Plot with Compensator');
% grid on;
% % saveas(gcf, 'figs/P3_BodePlot.png');
% 
% % Create a Nyquist Plot
% figure;
% hold on;
% nyquist(Lg);
% title('Loop Gain Nyquist Plot')
% grid on;
% 
% % Plot Poles and Zeros 
% figure;
% hold on;
% plot(real(CLZeros), imag(CLZeros), 'ro', 'MarkerSize', 8);
% plot(real(CLPoles), imag(CLPoles), 'rx', 'MarkerSize', 10);
% plot(real(OLZeros), imag(OLZeros), 'bo', 'MarkerSize', 8);
% plot(real(OLPoles), imag(OLPoles), 'bx', 'MarkerSize', 10);
% 
% xlabel('Real Axis');
% ylabel('Imaginary Axis');
% title('Pole-Zero Plot');
% grid on;
% axis equal; 
% legend('Closed Loop Zeros', 'Closed Loop Poles','Open Loop Zero','Open Loop Poles');
% hold off; 
% 
% 
% % Plot Tracking Error
% figure;
% semilogx(w, Ts_eval*100, 'LineWidth', 1.5);
% xlabel('Frequency (rad/s)');
% ylabel('Tracking Error (%)');
% title('Closed-Loop Tracking Error vs Frequency');
% grid on;
% % yline(5, '--r', '5% Error Threshold');
% ylim([0 100]);


%% Problem 4: Emperical Compensator
% Define Compensator
% Pole Cancellation: Cancel out the zero and pole with "approximate" cancellation
zeta_r = 0.2;
zeta_ar = 0.1;
omega_ar1 = 4.601;
omega_r1 = 8.22;
pole_cancel = (s^2 + 2*zeta_r*omega_r1*s + omega_r1^2)/ (s^2 + 2*zeta_ar*omega_ar1*s + omega_ar1^2);

% Proportional Gain
K = 1.2;

% Lead compensator
a1 = 0.4;
b1 = 1;
C1 = b1 / a1 * (s + a1) / (s + b1);

% Lag Compensator
a2 = 3;
b2 = 11;
C2 = b2 / a2 * (s + a2) / (s + b2);

% Total Compensator
C = K*C1*C2*pole_cancel;

% Evaluate the Compensator over Emperical Frequency Sweep 
C_eval = zeros(size(freq));
for i = 1:length(freq)
    C_eval(i) = evalfr(C, 1j * freq(i)); 
end

% Define Open Loop Magnitutde
OL_amp = exp(amp_db_P/20);
OL_mag = abs(OL_amp);

% Evalutate Loop Gain Magnitutde and Phase
Lg_Emperical = C_eval.*OL_mag;
Lg_mag = abs(Lg_Emperical);
Lg_phase = rad2deg(angle(Lg_Emperical));

% Evaluate Closed Loop Magnitutde and Phase
CL_Emperical = Lg_Emperical./(1+Lg_Emperical);
CL_mag = abs(CL_Emperical);
CL_phase = rad2deg(angle(CL_Emperical));

% Find Phase Margin
gain_crossing_idx = find(diff(sign(Lg_mag - 1)) ~= 0);
idx = gain_crossing_idx(1);
w1 = freq(idx);
w2 = freq(idx+1);
m1 = Lg_mag(idx);
m2 = Lg_mag(idx+1);
p1 = Lg_phase(idx);
p2 = Lg_phase(idx+1);
gain_crossover_freq = interp1([m1, m2], [w1, w2], 1);
phase_at_gc = interp1([w1, w2], [p1, p2], gain_crossover_freq);
PM = 180 + phase_at_gc;

% % Find the Gain Margin (RN it is INF)
% crossing_idx = find(Lg_phase <= -178); 
% idx = crossing_idx(1); 
% w1 = freq(idx-1);
% w2 = freq(idx+1);
% p1 = Lg_phase(idx-1);
% p2 = Lg_phase(idx+1);
% m1 = (Lg_mag(idx-1));
% m2 = (Lg_mag(idx+1));
% phase_crossover_freq = interp1([p1, p2], [w1, w2], -178);
% mag_at_pc = interp1([w1, w2], [m1, m2], phase_crossover_freq);
% GM = 20 * log10(mag_at_pc);

% Find the Closed Loop Bandwidth
idx = find(20*log10(CL_mag) < 20*log10(CL_mag(1))-3, 1);  
bw = freq(idx);

% Display Margins and CL Bandwidth
disp(['Gain Margin: ', num2str(GM), ' dB']);
disp(['Emperical Phase Margin: ', num2str(PM), ' deg']);
disp(['Empirical Closed-Loop Bandwidth: ', num2str(bw), ' rad/s']);

% Plot the Bode of the Loop Gain
figure;
subplot(2,1,1)
semilogx(freq, 20*log10(Lg_mag), 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('|-Lg(j\omega)| (dB)');
title('Magnitutde of the Empirical Loop-Gain');
grid on;
subplot(2,1,2)
semilogx(freq, Lg_phase, 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
title('Phase of Empirical Loop-Gain Response');
grid on;

% Plot the Bode of the Closed Loop
figure;
subplot(2,1,1)
semilogx(freq, 20*log10(CL_mag), 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('|T(j\omega)| (dB)');
title('Magnitutde of the Empirical Closed-Loop');
grid on;
subplot(2,1,2)
semilogx(freq, CL_phase, 'LineWidth', 1.5);
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
title('Phase of the Empirical Closed-Loop Response');
grid on;