% ASEN 5114 Project 1
% Authors: Carosn Kohlbenner, Thomas Dunnington, and Owen Craig4
clc;
clear;
close all;

fig_size = [100, 100, 800, 400];

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

% K = 25;
% params = {
%     %[2, 0.4, 0.015, omega_ar1, omega_ar1]; % anti-resonant notch
%     %[2, 0.035, 0.6, omega_r1, omega_r1]; % resonant notch
%     [2, 2*0.035, 0.5*0.015, omega_r1, omega_ar1]; % pole cancelation
%     [1, 1, 10]; % lead compensator
%     %[1, 1, 10]; % lead compensator
%     %[1, 0.3, 3]; % lag compensator
% };

% Parameters that describe the compensator
K = 44;
params = {
    [2, 5, 0.02, omega_ar1, omega_ar1]; % anti-resonant notch
    [2, 0.05, 1.1, omega_r1, omega_r1]; % resonant notch
};

% Calculate Gain and Phase Margin 
plot_compensators = true;
[Lg, C, fig] = loopgain(G, s, params, K, plot_compensators, "figs/Problem3/compensators.png");

% Plot the bode plot
[mag_g, phase_g, wout_g] = bode(G, freq_exp);
mag_g = squeeze(mag_g);
phase_g = squeeze(phase_g);

if plot_compensators
    figure(fig)
    fig.Position = fig_size;
    subplot(2,1,1)
    semilogx(wout_g, db(mag_g), 'k--', 'linewidth', 2, 'DisplayName', 'Empirical');
    hold on;
    title('Magnitude');
    xlabel('Frequency (rad/s)');
    ylabel('Amplitude (dB)');
    xlim([freq_exp(1), freq_exp(end)]);
    grid on;
    legend('show');

    % Plot the phase plot
    subplot(2,1,2)
    semilogx(wout_g, phase_g, 'k--', 'linewidth', 2, 'DisplayName', 'Empirical');
    hold on;
    title('Phase');
    xlabel('Frequency (rad/s)');
    ylabel('Phase (deg)');
    xlim([freq_exp(1), freq_exp(end)]);
    grid on;
    legend('show');
end

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

CL_bode_plot(Lg, "figs/Problem3/", "Prob3");

% Calculate the transfer function T = C/(1+Lg)
T_tf = C / (1 + Lg);

% Get Bode data for T_tf
[mag_T, phase_T, wout_T] = bode(T_tf);
mag_T = squeeze(mag_T);
phase_T = squeeze(phase_T);

% Create Bode plot for T_tf
figure_T_bode = figure;
figure_T_bode.Position = fig_size;

subplot(2,1,1);
semilogx(wout_T, db(mag_T));
title('Bode Plot Magnitude for T = C/(1+Lg)');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
grid on;

subplot(2,1,2);
semilogx(wout_T, phase_T);
title('Bode Plot Phase for T = C/(1+Lg)');
xlabel('Frequency (rad/s)');
ylabel('Phase (deg)');
grid on;
sgtitle('Bode Plot for T = C/(1+Lg)');


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

% figs = [figure, figure];

% a_min = 0.1;
% a_max = 2;
% for a = linspace(a_min, a_max, 10)

%     transparency = a/a_max;

%     mid_freq = 15;

%     K = 18;
%     params = {
%         [2, 5, 0.02, omega_ar1, omega_ar1]; % anti-resonant notch
%         [2, 0.05, 1.1, omega_r1, omega_r1]; % resonant notch
%         [1, 0.9, mid_freq]; % lead compensator
%         [1, 400, mid_freq]; % lag compensator
%     };

%     [~, C, ~] = loopgain(G, s, params, K, plot_compensators);

%     [c_mag, c_phase, c_wout] = bode(C, freq_exp);
%     c_mag = squeeze(c_mag);
%     c_phase = squeeze(c_phase);

%     Lg_exp.lg_mag = c_mag.*mag_exp;
%     Lg_exp.lg_phase = c_phase + rad2deg(phase_exp);
%     Lg_exp.lg_wout = c_wout;

%     Lg_exp.cl_mag = Lg_exp.lg_mag./(1+Lg_exp.lg_mag);
%     Lg_exp.cl_phase = Lg_exp.lg_phase - rad2deg(1./(1+Lg_exp.lg_mag));
%     Lg_exp.cl_wout = Lg_exp.lg_wout;

%     plot_lg(figs, Lg_exp, transparency, false);
% end

% MONEY %
K = 21;
params = {
    [2, 5, 0.02, omega_ar1, omega_ar1]; % anti-resonant notch
    [2, 0.05, 1.1, omega_r1, omega_r1]; % resonant notch
    [1, 1, 3]; % lead compensator
    [1, 400, 10]; % lag compensator
};

% Original
% K = 44;
% params = {
%     [2, 5, 0.02, omega_ar1, omega_ar1]; % anti-resonant notch
%     [2, 0.05, 1.1, omega_r1, omega_r1]; % resonant notch
% };

params_6 = params;
K6 = K;

% Calculate Gain and Phase Margin 
plot_compensators = true;
[~, C, fig] = loopgain(G, s, params, K, plot_compensators, "figs/Problem4/compensators.png");

if plot_compensators
    figure(fig)
    fig.Position = fig_size;
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
end

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
out_unsaturated = sim('model_step_unsaturated', 'StopTime', '100');  % Set simulation time to 100 seconds
% disp(pole(G))

step_plot = figure;
step_plot.Position = fig_size;
plot(out.input.Time, out.input.Data, 'r', 'LineWidth', 2, 'DisplayName', 'Reference Angle');
hold on;
plot(out.theta.Time, out.theta.Data, 'b', 'LineWidth', 2, 'DisplayName', 'With Saturation');
plot(out_unsaturated.theta.Time, out_unsaturated.theta.Data, 'b--', 'LineWidth', 1, 'DisplayName', 'Without Saturation');
title('Step Response');
xlabel('Time (s)');
ylabel('Angle (deg)');
grid on;
legend('show');
saveas(step_plot, 'figs/Problem5/step_response.png');

step_u_plot = figure;
step_u_plot.Position = fig_size;
plot(out.u.Time, out.u.Data, 'b', 'LineWidth', 2, 'DisplayName', 'With Saturation');
hold on;
plot(out_unsaturated.u.Time, out_unsaturated.u.Data, 'b--', 'LineWidth', 1, 'DisplayName', 'Without Saturation');
title('Step Response Controller Output');
xlabel('Time (s)');
ylabel('Torque (mNm)');
grid on;
legend('show');
saveas(step_u_plot, 'figs/Problem5/step_u_response.png');

% Sine response
w_in_set = [0.5, 2]; % rad/s
for i = 1:length(w_in_set)
    w_in = w_in_set(i);
    out(i) = sim('model_sine', 'StopTime', '100');  % Set simulation time to 100 seconds
    out_unsaturated(i) = sim('model_sine_unsaturated', 'StopTime', '100');  % Set simulation time to 100 seconds
    % disp(pole(G))
end

sine_plot = figure;
sine_plot.Position = fig_size;
for i = 1:length(w_in_set)
    subplot(length(w_in_set), 1, i);
    plot(out(i).input.Time, out(i).input.Data, 'r', 'LineWidth', 2, 'DisplayName', 'Reference Angle');
    hold on;
    plot(out(i).theta.Time, out(i).theta.Data, 'b', 'LineWidth', 2, 'DisplayName', 'With Saturation');
    plot(out_unsaturated(i).theta.Time, out_unsaturated(i).theta.Data, 'b--', 'LineWidth', 1, 'DisplayName', 'Without Saturation');
    title(sprintf('Sine Response - w_in = %g rad/s', w_in_set(i)));
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    grid on;
    legend('show');
end
saveas(sine_plot, 'figs/Problem5/sine_response.png');

u_plot = figure;
u_plot.Position = fig_size;
for i = 1:length(w_in_set)
    subplot(length(w_in_set), 1, i);
    plot(out(i).u.Time, out(i).u.Data, 'r', 'LineWidth', 2, 'DisplayName', 'With Saturation');
    hold on;
    plot(out_unsaturated(i).u.Time, out_unsaturated(i).u.Data, 'r--', 'LineWidth', 1, 'DisplayName', 'Without Saturation');
    title(sprintf('Sine Response - w_in = %g rad/s', w_in_set(i)));
    xlabel('Time (s)');
    ylabel('Torque (mNm)');
    grid on;
    legend('show');
end
saveas(u_plot, 'figs/Problem5/u_sine_response.png');

%%%%% Problem 6 %%%%%
disp('--- Problem 6 ---')
% prob6(C, s, "figs/Problem6/");
DC_gain = 10^(-15/20);
pole_1 = 0.3;
pole_2 = 0.4;
omega_ar1 = 4.58;   % anti-resonance frequency
omega_r1 = 2;    % resonance frequency
zeta_z = 0.015;
zeta_p = 0.035;
G_resonance = (s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2);
G = DC_gain * G_resonance * 1 / (s + pole_1) * 1 / (s + pole_2);

[Lg, C, fig] = loopgain(G, s, params_6, K6, true, "figs/Problem6/compensators.png");

CL_bode_plot(Lg, "figs/Problem6/", "Prob6");
[mag_g, phase_g, wout_g] = bode(G, freq_exp);
mag_g = squeeze(mag_g);
phase_g = squeeze(phase_g);

[mag_lg, phase_lg, wout_lg] = bode(Lg);
phase_lg = squeeze(phase_lg);
mag_lg = squeeze(mag_lg);

[mag_c, phase_c, wout_c] = bode(C);
phase_c = squeeze(phase_c);
mag_c = squeeze(mag_c);

if plot_compensators
    figure(fig)
    fig.Position = fig_size;
    subplot(2,1,1)
    semilogx(wout_g, db(mag_g), 'k--', 'linewidth', 2, 'DisplayName', 'Empirical');
    hold on;
    title('Magnitude');
    xlabel('Frequency (rad/s)');
    ylabel('Amplitude (dB)');
    xlim([freq_exp(1), freq_exp(end)]);
    grid on;
    legend('show');

    % Plot the phase plot
    subplot(2,1,2)
    semilogx(wout_g, phase_g, 'k--', 'linewidth', 2, 'DisplayName', 'Empirical');
    hold on;
    title('Phase');
    xlabel('Frequency (rad/s)');
    ylabel('Phase (deg)');
    xlim([freq_exp(1), freq_exp(end)]);
    grid on;
    legend('show');
end



function plot_lg(figs, Lg_exp, transparency, first_plot)

    mag_lg = Lg_exp.lg_mag;
    phase_lg = Lg_exp.lg_phase;
    wout_lg = Lg_exp.lg_wout;

    mag_cl = Lg_exp.cl_mag;
    phase_cl = Lg_exp.cl_phase;
    wout_cl = Lg_exp.cl_wout;

   
    % Calculate the closed loop bandwidth
    % Calculate the closed loop bandwidth for the experimental data
    closed_loop_bandwidth = wout_cl(find(db(mag_cl) <= -3, 1, 'first')); %Bandwidth is the frequency where the magnitude is -3 dB

    % Calculate the open loop stability margins
    % Phase margin
    gain_crossover_ind = find(diff(sign(db(mag_lg))) < 0 | diff(sign(db(mag_lg))) == 2);

    % Gain Margin
    phase_crossover_ind = find(diff(sign(phase_lg + 180)) < 0 | diff(sign(phase_lg + 180)) == 2);

    %% Loop Gain Bode Plot
    % Plot the model bode plot
    figure(figs(1));
    set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
    subplot(2,1,1)
    semilogx(wout_lg, db(mag_lg), 'b', 'linewidth', 2, 'color', [transparency, 0, 0]);
    hold on;
    if(~isempty(phase_crossover_ind))
        for i = 1:length(phase_crossover_ind)
            if abs(db(mag_lg(phase_crossover_ind(i)))) > 10
                line([wout_lg(phase_crossover_ind(i)), wout_lg(phase_crossover_ind(i))], [0, db(mag_lg(phase_crossover_ind(i)))], 'Color', 'g', 'LineWidth', 1.5);
            else
                line([wout_lg(phase_crossover_ind(i)), wout_lg(phase_crossover_ind(i))], [0, db(mag_lg(phase_crossover_ind(i)))], 'Color', 'r', 'LineWidth', 1.5);
            end
        end
    end
    yline(0, 'color', 'r', 'linestyle', ':', 'linewidth', 1.5)
    title('Magnitude');
    xlabel('Frequency (rad/s)');
    ylabel('Amplitude (dB)');
    xlim([wout_lg(1), wout_lg(end)]);
    grid on;

    % Plot the phase plot
    subplot(2,1,2)
    semilogx(wout_lg, phase_lg, 'linewidth', 2, 'color', [transparency, 0, 0]);
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
    figure(figs(2));
    set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
    subplot(2,1,1)
    semilogx(wout_cl, db(mag_cl), 'color', [transparency, 0, 0], 'linewidth', 2);
    hold on;
    xline(2*pi, 'k--', 'LineWidth', 1, 'DisplayName', '1 Hz', 'LabelVerticalAlignment', 'bottom'); % Add line at 1 Hz with label
    bandwidth_plot = xline(closed_loop_bandwidth, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--', 'label', closed_loop_bandwidth, 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'left');
    yline(-3, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
    title('Magnitude');
    xlabel('Frequency (rad/s)');
    ylabel('Amplitude (dB)');
    xlim([wout_cl(1), wout_cl(end)]);
    legend([bandwidth_plot], 'Bandwidth');
    grid on;

    % Plot the phase plot
    subplot(2,1,2)
    semilogx(wout_cl, phase_cl, 'color', [transparency, 0, 0], 'linewidth', 2);
    hold on;
    xline(2*pi, 'k--', 'LineWidth', 1, 'DisplayName', '1 Hz', 'Label', '1 Hz', 'LabelVerticalAlignment', 'bottom'); % Add line at 1 Hz with label
    % xline(2*pi, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
    title('Closed Loop Bode Plot');
    xlabel('Frequency (rad/s)');
    ylabel('Phase (deg)');
    xlim([wout_cl(1), wout_cl(end)]);
    sgtitle('Closed Loop Bode Plot')
    grid on;
end