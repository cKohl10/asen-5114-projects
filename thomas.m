% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig4
clc;
clear;
close all;

% Read the data from the file
data = readmatrix('data\Spacecraft_spin_module_frequency_response_data.xlsx');
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
semilogx(freq, amp_db_P, 'linewidth', 2);
hold on;
semilogx(wout, db(mag), 'r', 'linewidth', 2);
title('Bode Plot');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([freq(1), freq(end)]);
legend('Empirical', 'Analytical');
grid on;
saveas(gcf, 'figs/bode_plot.png');

% Plot the phase plot
figure;
semilogx(freq, rad2deg(phase_P), 'linewidth', 2);
hold on;
semilogx(freq, phase, 'color', 'r', 'linewidth', 2);
title('Phase Plot');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([freq(1), freq(end)]);
grid on;
legend('Empirical', 'Analytical');
saveas(gcf, 'figs/phase_plot.png');


%% Nyquist plotting
nyquist_exp(freq, amp_db_P, phase);



%% Problem 5
% Define controller transfer function
K = 1;
C = 0*s + K;

% Define the plant


% Step response
sinStepBool = 0; % Step
out = sim('model_sine_2024a', 'StopTime', '100');  % Set simulation time to 100 seconds

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
sinStepBool = 1;
w_in_set = [0.5, 2]; % rad/s
for i = 1:length(w_in_set)
    w_in = w_in_set(i);
    out(i) = sim('model_sine_2024a', 'StopTime', '100');  % Set simulation time to 100 seconds
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
