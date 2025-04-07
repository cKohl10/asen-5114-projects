% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig4
clc;
clear;
close all;

%% Plant
% Define the analytical transfer function
s = tf('s');
DC_gain = 10^(-15/20);
pole_1 = 0.3;
pole_2 = 0.4;
omega_ar1 = 4.58;   % anti-resonance frequency
omega_r1 = 8.347;    % resonance frequency
zeta_z = 0.015;
zeta_p = 0.035;
G_resonance = (s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2);
G = DC_gain * G_resonance * 1 / (s + pole_1) * 1 / (s + pole_2);

% Empirical Frequency Ranges
data = readmatrix('data\Spacecraft_spin_module_frequency_response_data.xlsx');
freq_exp = data(:,1)*2*pi;
mag_exp = data(:,2);
phase_exp = data(:,3);
w_min = min(freq_exp);
% w_min = 0;
w_max = max(freq_exp);

numG = G.Numerator{1};
denG = G.Denominator{1};
% [mag_PA, phase_PA, wout_PA] = bode(G, freq);
% mag_PA = squeeze(mag_PA);
% phase_PA = squeeze(phase_PA);
% % % Plot the bode plot
% bode_fig = figure;
% semilogx(freq, amp_db_PE);
% hold on;
% semilogx(wout_PA, db(mag_PA), 'r');
% title('Bode Plot');
% xlabel('Frequency (Hz)');
% ylabel('Amplitude (dB)');
% xlim([freq(1), freq(end)]);
% legend('P', 'P(s)');
% grid on;
% % saveas(gcf, 'figs/bode_plot.png');

% % Plot the phase plot
% phase_fig = figure;
% semilogx(freq, rad2deg(phase_PE));
% hold on;
% semilogx(wout_PA, phase_PA, 'r');
% title('Phase Plot');
% xlabel('Frequency (Hz)');
% ylabel('Phase (deg)');
% xlim([freq(1), freq(end)]);
% grid on;
% legend('P', 'P(s)');
% % saveas(gcf, 'figs/phase_plot.png');

%%%%% Problem 3 %%%%%
disp('--- Problem 3 ---')

% type = 1: Lag Compensator / Lead Compensator
% type = 2: Notch Filter
% type = 3: Proportional Gain
% type = 4: Integral Gain
% type = 5: Derivative Gain

K = 25;
params = {
    %[2, 0.4, 0.015, omega_ar1, omega_ar1]; % anti-resonant notch
    %[2, 0.035, 0.6, omega_r1, omega_r1]; % resonant notch
    [2, 2*0.035, 0.5*0.015, omega_r1, omega_ar1]; % pole cancelation
    [1, 1, 10]; % lead compensator
    %[1, 1, 10]; % lead compensator
    %[1, 0.3, 3]; % lag compensator
};

% Calculate Gain and Phase Margin 
plot_compensators = false;
[Lg, C] = loopgain(G, s, params, K, plot_compensators);

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
disp(['Gain Margin: ', num2str(GM), 'dB']);
disp(['Phase Margin: ', num2str(PM), 'deg']);
disp(['Closed-loop bandwidth: ', num2str(bw), ' rad/s']);

CL_bode_plot(Lg, "figs/Problem3/", "Prob3");

numCA = C.Numerator{1};
denCA = C.Denominator{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Problem 4 %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--- Problem 4 ---')

% type = 1: Lag Compensator / Lead Compensator
% type = 2: Notch Filter
% type = 3: Proportional Gain
% type = 4: Integral Gain
% type = 5: Derivative Gain

K = 35;
params = {
    [2, 5, 0.02, omega_ar1, omega_ar1]; % anti-resonant notch
    [2, 0.1, 1.1, omega_r1, omega_r1]; % resonant notch
    % [2, 2*0.035, 0.5*0.015, omega_r1, omega_ar1]; % pole cancelation
    [1, 7, 10]; % lead compensator
    [1, 1, 5]; % lead compensator
    %[1, 8, 6]; % lag compensator
    [1, 40, 2]; % lag compensator
    % [6, 2*pi]; % low pass filter
};

% Calculate Gain and Phase Margin 
plot_compensators = true;
[~, C, fig] = loopgain(G, s, params, K, plot_compensators);

figure(fig)
fig.Position = [100, 100, 1000, 800];
subplot(2,1,1)
semilogx(freq_exp, db(mag_exp), 'k--', 'linewidth', 2, 'DisplayName', 'Empirical');
hold on;
title('Magnitude');
xlabel('Frequency (rad/s)');
ylabel('Amplitude (dB)');
xlim([freq_exp(1), freq_exp(end)]);
grid on;
legend('show');

% Plot the phase plot
subplot(2,1,2)
semilogx(freq_exp, rad2deg(phase_exp), 'k--', 'linewidth', 2, 'DisplayName', 'Empirical');
hold on;
title('Phase');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
xlim([freq_exp(1), freq_exp(end)]);
grid on;
legend('show');

[c_mag, c_phase, c_wout] = bode(C, freq_exp);
c_mag = squeeze(c_mag);
c_phase = squeeze(c_phase);

Lg_exp.lg_mag = c_mag.*mag_exp;
Lg_exp.lg_phase = c_phase + rad2deg(phase_exp);
Lg_exp.lg_wout = c_wout;

Lg_exp.cl_mag = Lg_exp.lg_mag./(1+Lg_exp.lg_mag);
Lg_exp.cl_phase = Lg_exp.lg_phase - rad2deg(1./(1+Lg_exp.lg_mag));
Lg_exp.cl_wout = Lg_exp.lg_wout;

CL_bode_plot(Lg_exp, "figs/Problem4/", "Prob4");

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
saveas(gcf, 'figs/Problem5/step_response.png');

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
saveas(gcf, 'figs/Problem5/sine_response.png');

figure;
for i = 1:length(w_in_set)
    subplot(length(w_in_set), 1, i);
    plot(out(i).u.Time, out(i).u.Data, 'r--');
    hold on;
    title(sprintf('Sine Response - w_in = %g rad/s', w_in_set(i)));
    xlabel('Time (s)');
    ylabel('Torque input (Nm)');
    grid on;
end
saveas(gcf, 'figs/Problem5/u_sine_response.png');

%%%%% Problem 6 %%%%%
% prob6(C, s)
