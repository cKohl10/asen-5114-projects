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
numG = G.Numerator{1};
denG = G.Denominator{1};
[mag, phase, wout] = bode(G, freq);
mag = squeeze(mag);

% Plot the bode plot
figure;
semilogx(freq, amp_db_P);
hold on;
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
title('Phase Plot');
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
xlim([freq(1), freq(end)]);
grid on;
legend('P', 'P_0');
saveas(gcf, 'figs/phase_plot.png');


%%%%% Problem 2 %%%%%

% % Simple controller - plot all K values on one Nyquist plot
% K = [0.1, 1, 10, 100, 1000];
% figure;
% hold on;
% legendEntries = cell(1, length(K));

% for i = 1:length(K)
%     C = K(i);
%     nyquist(C*G);
%     legendEntries{i} = sprintf('K = %g', K(i));
% end

% title('Nyquist Plot - Multiple Gain Values');
% xlabel('Real Axis');
% ylabel('Imaginary Axis');
% legend(legendEntries);
% grid on;
% saveas(gcf, 'figs/nyquist_plot_combined.png');

%%%%% Problem 3 %%%%%
cz1 = 1;
cp1 = 2;
C1 = (s+cz1)/(s+cp1);

K_c = 1;
C = C1*K_c;

numC = C.Numerator{1};
denC = C.Denominator{1};

%%%%% Problem 5 %%%%%

% Step response
out = sim('model_step', 'StopTime', '100');  % Set simulation time to 100 seconds
disp(pole(G))

figure;
plot(out.input.Time, out.input.Data, 'r--');
hold on;
plot(out.theta.Time, out.theta.Data);
title('Step Response');
xlabel('Time (s)');
ylabel('Theta');
grid on;
saveas(gcf, 'figs/step_response.png');

% Sine response
w_in_set = [0.5, 2]; % rad/s
for i = 1:length(w_in_set)
    w_in = w_in_set(i);
    out(i) = sim('model_sine', 'StopTime', '100');  % Set simulation time to 100 seconds
    disp(pole(G))
end

figure;
for i = 1:length(w_in_set)
    subplot(length(w_in_set), 1, i);
    plot(out(i).input.Time, out(i).input.Data, 'r--');
    hold on;
    plot(out(i).theta.Time, out(i).theta.Data);
    title(sprintf('Sine Response - w_in = %g rad/s', w_in_set(i)));
    xlabel('Time (s)');
    ylabel('Theta');
    grid on;
end
saveas(gcf, 'figs/sine_response.png');
