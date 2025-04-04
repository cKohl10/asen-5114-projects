% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig4
clc;
clear;
close all;

% Read the data from the file
data = table2array(readtable('data\Spacecraft_spin_module_frequency_response_data.xlsx'));
freq = data(:,1);
mag_PE = data(:,2);
amp_db_PE = 20*log10(mag_PE) - 20*log10(freq);
amp_db_P0 = 20*log10(mag_PE);
phase_PE = data(:,3) - pi/2;
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
[mag_PA, phase_PA, wout_PA] = bode(G, freq);
mag_PA = squeeze(mag_PA);
phase_PA = squeeze(phase_PA);
% Plot the bode plot
bode_fig = figure;
semilogx(freq, amp_db_PE);
hold on;
semilogx(wout_PA, db(mag_PA), 'r');
title('Bode Plot');
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
xlim([freq(1), freq(end)]);
legend('P', 'P(s)');
grid on;
saveas(gcf, 'figs/bode_plot.png');

% Plot the phase plot
phase_fig = figure;
semilogx(freq, rad2deg(phase_PE));
hold on;
semilogx(wout_PA, phase_PA, 'r');
title('Phase Plot');
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
xlim([freq(1), freq(end)]);
grid on;
legend('P', 'P(s)');
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
disp('--- Problem 3 ---')

% Lead Compensator
% cz1 = 0.1;
% cp1 = 0.5;
% C1 = (s/cz1 + 1)/(s/cp1 + 1);
% K_c = 100;

% % Pole Cancelation resonance pair
% c_omega_ar = omega_r1;
% c_omega_r = omega_ar1;
% C1 = (s^2 + 2*zeta_z*c_omega_ar*s + c_omega_ar^2) / (s^2 + 2*zeta_p*c_omega_r*s + c_omega_r^2);
% K_c1 = 1;

% % Pole Cancelation slope
% cz1 = 0.1;
% cp1 = 1;
% C2 = (s/cz1 + 1)/(s/cp1 + 1);
% K_c2 = 40;

% C = C1*K_c1*C2*K_c2;

% Define Compensator Properties
z1 = 0.01;   
p1 = 6*pi; 
z2 = 1;
p2 = 6*pi;

% Compensator Gain (Tuned for Stability)
K = 80000;
zeta_z = 0.01;
zeta_p = 0.06;

% Combined Notch Filter to Damp Resonance and Anti-Resonance
notch = (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2)/(s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2);
C = K * (s + z1)/(s + p1) * (s + z2)/(s + p2) * 1/s * notch;
% Calculate Gain and Phase Margin 
Lg = -C*G;

[c_mag, c_phase, c_wout] = bode(C, freq);
[cg_mag, cg_phase, cg_wout] = bode(C*G, freq);
cg_mag = squeeze(cg_mag);
cg_phase = squeeze(cg_phase);
c_mag = squeeze(c_mag);
c_phase = squeeze(c_phase);

% Find zero crossings on the bode plot with linear interpolation
zero_crossings = [];
for i = 1:length(cg_mag)-1
    if (cg_mag(i) - 1) * (cg_mag(i+1) - 1) <= 0 % Check if there is a sign change
        % Linear interpolation between points
        x1 = freq(i);
        x2 = freq(i+1);
        y1 = cg_mag(i) - 1;
        y2 = cg_mag(i+1) - 1;
        
        % Calculate exact crossing point
        crossover_freq = x1 - y1*(x2 - x1)/(y2 - y1);
        zero_crossings = [zero_crossings, crossover_freq];
    end
end

disp('Zero crossing frequencies (Hz):');
disp(zero_crossings);

% Get the phase margin and gain margin at the zero crossings
if ~isempty(zero_crossings)
    for i = 1:length(zero_crossings)
        % Get the closest phase value in the fitted plot to the zero crossing frequency
        [~, idx] = min(abs(freq - zero_crossings(i)));
        zc_phase(i) = cg_phase(idx);
        pm(i) = zc_phase(i) + 180;
    end
    disp('Phase margin at zero crossings (deg):');
    disp(pm);
else 
    disp('No zero crossings found');
    pm = NaN;
end

% Find -180 degree crossings with linear interpolation
phase_180_crossings = [];
for i = 1:length(cg_phase)-1
    if (cg_phase(i) + 180) * (cg_phase(i+1) + 180) <= 0 % Check if there is a sign change around -180
        % Linear interpolation between points
        x1 = freq(i);
        x2 = freq(i+1);
        y1 = cg_phase(i) + 180;
        y2 = cg_phase(i+1) + 180;
        
        % Calculate exact crossing point
        crossover_freq = x1 - y1*(x2 - x1)/(y2 - y1);
        phase_180_crossings = [phase_180_crossings, crossover_freq];
    end
end

if isempty(phase_180_crossings)
    disp('Infinite gain margin');
    gm = NaN;
else
    % Calculate gain margins at each -180 degree crossing
    for i = 1:length(phase_180_crossings)
        [~, idx] = min(abs(freq - phase_180_crossings(i)));
        gm(i) = -db(cg_mag(idx));  % Gain margin in dB
    end
    disp('Gain margins (dB):');
    disp(gm);
end

% Find -3dB crossing (bandwidth)
bandwidth_idx = find(db(cg_mag) <= -3, 1, 'first');
if ~isempty(bandwidth_idx)
    bandwidth_freq = freq(bandwidth_idx);
    disp(['Bandwidth = ', num2str(bandwidth_freq), ' Hz']);
end

% Plot the bode plot
f_bode_3 = figure();
f_bode_3.Position = [100, 100, 800, 600];
semilogx(wout_PA, db(mag_PA), 'r');
hold on;
semilogx(c_wout, db(c_mag), 'b--');
semilogx(cg_wout, db(cg_mag), 'g');
if ~isempty(phase_180_crossings)
    scatter(phase_180_crossings, -gm, 'r', 'o');
end
yline(0, 'k--');
for i = 1:length(zero_crossings)
    xline(zero_crossings(i), 'k--');
end
grid on;
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
title('Problem 3 Bode Plot');
legend('P(s)', 'C(s)', 'C(s)*P(s)', 'Gain Margin Points', 'Location', 'eastoutside');
saveas(gcf, 'figs/bode_plot.png');

% Plot the phase plot
f_phase_3 = figure(phase_fig);
f_phase_3.Position = f_bode_3.Position;
hold on;
semilogx(wout_PA, phase_PA, 'r');
semilogx(c_wout, c_phase, 'b--');
semilogx(cg_wout, cg_phase, 'g');
if ~isempty(zero_crossings)
    scatter(zero_crossings, zc_phase, 'r', 'o');
end
yline(-180, 'k--');
ylabel('Phase (deg)');
xlabel('Frequency (Hz)');
title('Problem 3 Phase Plot');
legend('P(s)', 'C(s)', 'C(s)*P(s)', 'Phase Margin Points', 'Location', 'eastoutside');
saveas(gcf, 'figs/phase_plot.png');

numCA = C.Numerator{1};
denCA = C.Denominator{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Problem 4 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--- Problem 4 ---')

% Lead Compensator
% cz1 = 0.1;
% cp1 = 0.5;
% C1 = (s/cz1 + 1)/(s/cp1 + 1);
% K_c = 100;

% Pole Cancelation resonance pair
c_omega_ar = omega_r1;
c_omega_r = omega_ar1;
C1 = (s^2 + 2*zeta_z*c_omega_ar*s + c_omega_ar^2) / (s^2 + 2*zeta_p*c_omega_r*s + c_omega_r^2);
K_c1 = 1;

% Pole Cancelation slope
cz1 = 0.1;
cp1 = 1;
C2 = (s/cz1 + 1)/(s/cp1 + 1);
K_c2 = 40;

C = C1*K_c1*C2*K_c2;
[c_mag, c_phase, c_wout] = bode(C, freq);
c_mag = squeeze(c_mag);
c_phase = squeeze(c_phase);

cg_mag = c_mag.*mag_PE;
cg_phase = c_phase + rad2deg(phase_PE);

% Find zero crossings on the bode plot with linear interpolation
zero_crossings = [];
for i = 1:length(cg_mag)-1
    if (cg_mag(i) - 1) * (cg_mag(i+1) - 1) <= 0 % Check if there is a sign change
        % Linear interpolation between points
        x1 = freq(i);
        x2 = freq(i+1);
        y1 = cg_mag(i) - 1;
        y2 = cg_mag(i+1) - 1;
        
        % Calculate exact crossing point
        crossover_freq = x1 - y1*(x2 - x1)/(y2 - y1);
        zero_crossings = [zero_crossings, crossover_freq];
    end
end

disp('Zero crossing frequencies (Hz):');
disp(zero_crossings);

% Get the phase margin and gain margin at the zero crossings
if ~isempty(zero_crossings)
    for i = 1:length(zero_crossings)
        % Get the closest phase value in the fitted plot to the zero crossing frequency
        [~, idx] = min(abs(freq - zero_crossings(i)));
        zc_phase(i) = cg_phase(idx);
        pm(i) = zc_phase(i) + 180;
    end
    disp('Phase margin at zero crossings (deg):');
    disp(pm);
else 
    disp('No zero crossings found');
    pm = NaN;
end

% Find -180 degree crossings with linear interpolation
phase_180_crossings = [];
for i = 1:length(cg_phase)-1
    if (cg_phase(i) + 180) * (cg_phase(i+1) + 180) <= 0 % Check if there is a sign change around -180
        % Linear interpolation between points
        x1 = freq(i);
        x2 = freq(i+1);
        y1 = cg_phase(i) + 180;
        y2 = cg_phase(i+1) + 180;
        
        % Calculate exact crossing point
        crossover_freq = x1 - y1*(x2 - x1)/(y2 - y1);
        phase_180_crossings = [phase_180_crossings, crossover_freq];
    end
end

if isempty(phase_180_crossings)
    disp('Infinite gain margin');
    gm = NaN;
else
    % Calculate gain margins at each -180 degree crossing
    for i = 1:length(phase_180_crossings)
        [~, idx] = min(abs(freq - phase_180_crossings(i)));
        gm(i) = -db(cg_mag(idx));  % Gain margin in dB
    end
    disp('Gain margins (dB):');
    disp(gm);
end

% Find -3dB crossing (bandwidth)
bandwidth_idx = find(db(cg_mag) <= -3, 1, 'first');
if ~isempty(bandwidth_idx)
    bandwidth_freq = freq(bandwidth_idx);
    disp(['Bandwidth = ', num2str(bandwidth_freq), ' Hz']);
end

% Plot the bode plot
f_bode_4 = figure();
f_bode_4.Position = [100, 100, 800, 600];
semilogx(freq, db(mag_PE), 'r');
hold on;
semilogx(c_wout, db(c_mag), 'b--');
semilogx(cg_wout, db(cg_mag), 'g');
if ~isempty(phase_180_crossings)
    scatter(phase_180_crossings, -gm, 'r', 'o');
end
yline(0, 'k--');
for i = 1:length(zero_crossings)
    xline(zero_crossings(i), 'k--');
end
grid on;
ylabel('Magnitude (dB)');
xlabel('Frequency (Hz)');
title('Problem 4 Bode Plot');
legend('P(s)', 'C(s)', 'C(s)*P(s)', 'Gain Margin Points', 'Location', 'eastoutside');
saveas(gcf, 'figs/bode_plot.png');

% Plot the phase plot
f_phase_4 = figure();
f_phase_4.Position = f_bode_4.Position;
semilogx(freq, rad2deg(phase_PE), 'r');
hold on;
semilogx(c_wout, c_phase, 'b--');
semilogx(cg_wout, cg_phase, 'g');
% if ~isempty(zero_crossings)
%     scatter(zero_crossings, zc_phase, 'r', 'o');
% end
yline(-180, 'k--');
ylabel('Phase (deg)');
xlabel('Frequency (Hz)');
title('Problem 4 Phase Plot');
legend('P(s)', 'C(s)', 'C(s)*P(s)', 'Phase Margin Points', 'Location', 'eastoutside');
grid on;
saveas(gcf, 'figs/phase_plot.png');

numCE = C.Numerator{1};
denCE = C.Denominator{1};


%%%%% Problem 5 %%%%%

numC = numCA;
denC = denCA;

% Step response
out = sim('model_step', 'StopTime', '100');  % Set simulation time to 100 seconds
% disp(pole(G))

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
    % disp(pole(G))
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

figure;
for i = 1:length(w_in_set)
    subplot(length(w_in_set), 1, i);
    plot(out(i).u.Time, out(i).u.Data, 'r--');
    hold on;
    title(sprintf('Sine Response - w_in = %g rad/s', w_in_set(i)));
    xlabel('Time (s)');
    ylabel('u');
    grid on;
end
% saveas(gcf, 'figs/sine_response.png');
