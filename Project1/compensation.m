%% Compensation Design
% This script designs a compensator to achieve design objectives
close all; clear; clc;

%% Plant
% Define the analytical transfer function
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


% Empirical Frequency Ranges
data = readmatrix('data/Spacecraft_spin_module_frequency_response_data.xlsx');
freq_exp = data(:,1)*2*pi;
w_min = min(freq_exp);
% w_min = 0;
w_max = max(freq_exp);

%% Compensation
%%%% Anti Notch Filter to Cancel anti-resonance
w_n = omega_ar1;
w_d = omega_ar1;
zeta_n = 0.4;     
zeta_d = 0.015; 

% Notch Filter Transfer Function
notch_filter_ar = (s^2 + 2*zeta_n*w_n*s + w_n^2) / (s^2 + 2*zeta_d*w_d*s + w_d^2);


%%%% Notch Filter to cancel out resonance
w_n = omega_r1;
w_d = omega_r1;
zeta_n = 0.6;     
zeta_d = 0.035; 

% Notch Filter Transfer Function
notch_filter_r = (s^2 + 2*zeta_d*w_d*s + w_d^2) / (s^2 + 2*zeta_n*w_n*s + w_n^2);

%%%% Pole Cancellation
% Cancel out the zero and pole with "approximate" cancellation
zeta_r = 2*0.035;
zeta_ar = 0.5*0.015;
pole_cancel = (s^2 + 2*zeta_r*omega_r1*s + omega_r1^2) / (s^2 + 2*zeta_ar*omega_ar1*s + omega_ar1^2);

% Visulize the cancellation
[mag_n, phase_n, wout_n] = bode(pole_cancel);
phase_n = squeeze(phase_n);
mag_n = squeeze(mag_n);

% Proportional Gain
K = 25;

%%%% Lead compensator
a1 = 1;
b1 = 10;
C1 = b1 / a1 * (s + a1) / (s + b1);

% Bode of lead
[mag_lead, phase_lead, wout_lead] = bode(C1);
phase_lead = squeeze(phase_lead);
mag_lead = squeeze(mag_lead);

% Second Compensator
a2 = 1;
b2 = 10;
C2 = b2 / a2 * (s + a2) / (s + b2);

% Third Compensator (Lag)
% a3 = 20;
% b3 = 4;
% C3 = b3 / a3 * (s + a3) / (s + b3);
a3 = 0.3;
b3 = 3;
C3 = b3 / a3 * (s + a3) / (s + b3);

%Controller
% C1 = 1;
C2 = 1;
C3 = 1;
notch_filter_ar = 1;
notch_filter_r = 1;
%pole_cancel = 1;
% K = 1;
C = K*C1*C2*C3*notch_filter_ar*notch_filter_r*pole_cancel;

% Negative Loop Gain
Lg = C*G;

% Compute the frequency response
[mag_lg, phase_lg, wout_lg] = bode(Lg);
phase_lg = squeeze(phase_lg);
mag_lg = squeeze(mag_lg);

% Calculate the closed loop response
T = feedback(Lg, 1);
[mag_cl, phase_cl, wout_cl] = bode(T);
phase_cl = squeeze(phase_cl);
mag_cl = squeeze(mag_cl);

% Calculate the closed loop bandwidth
closed_loop_bandwidth = bandwidth(T);

% Calculate the open loop stability margins
% Phase margin
gain_crossover_ind = find(diff(sign(db(mag_lg))) < 0 | diff(sign(db(mag_lg))) == 2);

% Gain Margin
phase_crossover_ind = find(diff(sign(phase_lg + 180)) < 0 | diff(sign(phase_lg + 180)) == 2);

%% Loop Gain Bode Plot
% Plot the model bode plot
figure;
set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
subplot(2,1,1)
semilogx(wout_lg, db(mag_lg), 'b', 'linewidth', 2);
hold on;
yline(-10, 'color', 'g', 'linestyle', ':', 'linewidth', 1.5)
if(~isempty(phase_crossover_ind))
    for i = 1:length(phase_crossover_ind)
        if db(mag_lg(phase_crossover_ind(i))) < -10
            good_gain = line([wout_lg(phase_crossover_ind(i)), wout_lg(phase_crossover_ind(i))], [0, db(mag_lg(phase_crossover_ind(i)))], 'Color', 'g', 'LineWidth', 1.5);
        else
            bad_gain = line([wout_lg(phase_crossover_ind(i)), wout_lg(phase_crossover_ind(i))], [0, db(mag_lg(phase_crossover_ind(i)))], 'Color', 'r', 'LineWidth', 1.5);
        end
    end
end
yline(0, 'color', 'r', 'linestyle', ':', 'linewidth', 1.5)
title('Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([wout_lg(1), wout_lg(end)]);
grid on;

% Add legend
try
    legend([good_gain, bad_gain], 'Good Gain Margin', 'Bad Gain Margin', 'location' ,'best');
catch
    try 
        legend([good_gain], 'Good Gain Margin', 'location' ,'best');
    catch
        try
            legend([bad_gain], 'Bad Gain Margin', 'location' ,'best');
        catch
            x = 1;
        end
    end
end

% Plot the phase plot
subplot(2,1,2)
semilogx(wout_lg, phase_lg, 'linewidth', 2);
hold on;
for i = 1:length(gain_crossover_ind)
    % Plot green if good phase margin, otherwise red
    if phase_lg(gain_crossover_ind(i)) >= -140
        good_phase = line([wout_lg(gain_crossover_ind(i)), wout_lg(gain_crossover_ind(i))], [-180, phase_lg(gain_crossover_ind(i))], 'Color', 'g', 'LineWidth', 1.5);
    else
        bad_phase = line([wout_lg(gain_crossover_ind(i)), wout_lg(gain_crossover_ind(i))], [-180, phase_lg(gain_crossover_ind(i))], 'Color', 'r', 'LineWidth', 1.5);
    end
end
yline(-180, 'color', 'r', 'linestyle', ':', 'linewidth', 1.5)
yline(-140, 'color', 'g', 'linestyle', ':', 'linewidth', 1.5)
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([wout_lg(1), wout_lg(end)]);
sgtitle('Loop Gain Bode Plot')
grid on;

% Add legend
try
    legend([good_phase, bad_phase], 'Good Phase Margin', 'Bad Phase Margin', 'location' ,'best');
catch
    try 
        legend([good_phase], 'Good Phase Margin', 'location' ,'best');
    catch
        try
            legend([bad_phase], 'Bad Phase Margin', 'location' ,'best');
        catch
            x = 1;
        end
    end
end


%% Closed Loop Bode Plot
figure;
set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
subplot(2,1,1)
semilogx(wout_cl, db(mag_cl), 'color', 'g', 'linewidth', 2);
hold on;
bandwidth_plot = xline(closed_loop_bandwidth, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--', 'label', closed_loop_bandwidth, 'LabelVerticalAlignment', 'bottom');
yline(-3, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
title('Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([wout_cl(1), wout_cl(end)]);
legend([bandwidth_plot], 'Bandwidth');
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(wout_cl, phase_cl, 'color', 'g', 'linewidth', 2);
hold on;
% xline(2*pi, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([wout_cl(1), wout_cl(end)]);
sgtitle('Closed Loop Bode Plot')
grid on;


%% Nyquist Plot
% Create a nyquist plot of the model
figure();
nyquist(Lg, {w_min, w_max});
axis equal;


%% Notch Filter Plot
figure;
set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
subplot(2,1,1)
semilogx(wout_n, db(mag_n), 'b', 'linewidth', 2);
hold on;
title('Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([wout_n(1), wout_n(end)]);
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(wout_n, phase_n, 'b', 'linewidth', 2);
hold on;
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([wout_n(1), wout_n(end)]);
sgtitle('Pole Zero Cancellation Compensator Bode Plot')
grid on;


%% Lead Plot
% Just the lead compensator Bode plot
figure;
set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
subplot(2,1,1)
semilogx(wout_lead, db(mag_lead), 'b', 'linewidth', 2);
hold on;
title('Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([wout_lead(1), wout_lead(end)]);
grid on;

% Plot the phase plot
subplot(2,1,2)
semilogx(wout_lead, phase_lead, 'b', 'linewidth', 2);
hold on;
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([wout_lead(1), wout_lead(end)]);
sgtitle('Lead Compensator Bode Plot')
grid on;

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
    Lg_eval(i) = evalfr(Lg, w(i));
    Ts_eval(i) = 1/abs(1-Lg_eval(i)); 
end

% Display the margins and bandwidth 
disp(' --- Problem 3 --- ')
disp(['Gain Margin: ', num2str(GM), 'dB']);
disp(['Phase Margin: ', num2str(PM), 'deg']);
disp(['Closed-loop bandwidth: ', num2str(bw), ' rad/s']);