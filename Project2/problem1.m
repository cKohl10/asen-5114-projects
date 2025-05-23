%% Main Script for Problem 1
% Authors: Thomas Dunnington, Owen Craig, Carson Kohbrenner
close all; clear; clc;


%% Frequency Response System Identification
% Read the data from the file
data = readmatrix('data/Spacecraft_spin_module_frequency_response_data.xlsx');
freq_exp = data(:,1)*2*pi;
amp = data(:,2);
amp_db_P = 20*log10(amp) - 20*log10(freq_exp);
phase_P = data(:,3) - pi/2;


% Transfer function with anti resonance, resonance pair
s = tf('s');
DC_gain = 10^(-15/20);
pole_1 = 0;
pole_2 = 0.9;
omega_ar1 = 4.601;   % anti-resonance frequency
omega_r1 = 8.347;    % resonance frequency
zeta_z = 0.015;
zeta_p = 0.035;
G_resonance = (s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2);
G = DC_gain * G_resonance * 1 / (s + pole_1) * 1 / (s + pole_2);

% Compute the frequency response
[mag, phase, wout] = bode(G, {min(freq_exp), max(freq_exp)});
phase = squeeze(phase);
mag = squeeze(mag);

%% Try Using tfest
% Truncate to below 3 Hz
amp_db_P_tfest = amp_db_P(freq_exp <= 3*2*pi);
phase_P_tfest = phase_P(freq_exp <= 3*2*pi);
freq_exp_tfest = freq_exp(freq_exp <= 3*2*pi);

% Convert to complex response
mag_tfest = 10.^(amp_db_P_tfest ./ 20);
resp = mag_tfest .* exp(1j * phase_P_tfest);

% Group and average data by unique frequencies
[unique_freqs, ~, idx] = unique(freq_exp_tfest);

% Preallocate complex response
avg_resp = zeros(size(unique_freqs));

for i = 1:length(unique_freqs)
    avg_resp(i) = mean(resp(idx == i));
end

% Create frequency response data object
Ts = 0;  % Continuous-time system
frd_data = idfrd(avg_resp(:), unique_freqs(:), Ts);

% Estimate the transfer function
np = 4; % number of poles
nz = 2; % number of zeros
tfest_sys = tfest(frd_data, np, nz);

% Compute the frequency response
[mag_tfest, phase_tfest, wout_tfest] = bode(tfest_sys, {min(freq_exp), max(freq_exp)});
phase_tfest = squeeze(phase_tfest);
mag_tfest = squeeze(mag_tfest);


%% Create Minimal State Space Representation
% Controllable canonical form
[A, B, C, D] = tf2ss(G.Numerator{1}, G.Denominator{1});

% Get observable canonical form
A_obs = A';
B_obs = C';
C_obs = B';
D_obs = D';

% Save state space realizations
save('Data/state_space.mat' , 'A', 'B', 'C', 'D', 'A_obs', 'B_obs', 'C_obs', 'D_obs');

% Check controllability
ctrl_mat = ctrb(A, B);

% Check observability
obsv_mat = obsv(A, C);


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
xline(3*2*pi, 'linestyle', '--', 'color', 'r', 'label', 'Model Cutoff')
title('Magnitude');
ylabel('Amplitude (dB)');
xlim([freq_exp(1), freq_exp(end)]);
legend('Empirical', 'Analytical', 'location' , 'best');
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(freq_exp, rad2deg(phase_P), 'linewidth', 2);
hold on;
semilogx(wout, phase, 'color', 'r', 'linewidth', 2);
xline(3*2*pi, 'linestyle', '--', 'color', 'r', 'label', 'Model Cutoff')
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([freq_exp(1), freq_exp(end)]);
sgtitle('Bode Plot')
grid on;
