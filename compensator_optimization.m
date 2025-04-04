%% Compensator Optimization
% This script optimizes a compensator to meet the following criteria:
% 1. 1 Hz bandwidth (2π rad/s)
% 2. At least 40 degree phase margin
% 3. At least 10 dB gain margin

close all; clear; clc;

%% Plant Definition
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

% Target bandwidth in rad/s (1 Hz = 2π rad/s)
target_bandwidth = 2*pi;

%% Optimization Setup
% Define the parameters to optimize
% x = [K, a1, b1, a2, b2, a3, b3, zeta_n_ar, zeta_d_ar, zeta_n_r, zeta_d_r]

% Initial values
x0 = [150, 7, 13, 2, 20, 10, 20, 0.4, 0.015, 0.6, 0.035];

% Lower and upper bounds
lb = [1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01, 0.01];
ub = [1000, 30, 50, 30, 50, 30, 50, 1, 1, 1, 1];

% Optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 1000);

% Run the optimization
[x_opt, fval] = fmincon(@(x) compensator_cost(x, G, target_bandwidth), x0, [], [], [], [], lb, ub, [], options);

% Display results
disp('Optimized Parameters:');
disp(['K = ', num2str(x_opt(1))]);
disp(['a1 = ', num2str(x_opt(2)), ', b1 = ', num2str(x_opt(3))]);
disp(['a2 = ', num2str(x_opt(4)), ', b2 = ', num2str(x_opt(5))]);
disp(['a3 = ', num2str(x_opt(6)), ', b3 = ', num2str(x_opt(7))]);
disp(['zeta_n_ar = ', num2str(x_opt(8)), ', zeta_d_ar = ', num2str(x_opt(9))]);
disp(['zeta_n_r = ', num2str(x_opt(10)), ', zeta_d_r = ', num2str(x_opt(11))]);

% Evaluate and display performance metrics with optimized parameters
[cost, metrics] = compensator_cost(x_opt, G, target_bandwidth, true);
disp('Performance Metrics:');
disp(['Bandwidth = ', num2str(metrics.bandwidth), ' rad/s (', num2str(metrics.bandwidth/(2*pi)), ' Hz)']);
disp(['Phase Margin = ', num2str(metrics.phase_margin), ' degrees']);
disp(['Gain Margin = ', num2str(metrics.gain_margin), ' dB']);

% Plot the optimized system
plot_system(x_opt, G, target_bandwidth);

%% Cost Function
function [cost, metrics] = compensator_cost(x, G, target_bandwidth, return_metrics)
    % Default for return_metrics
    if nargin < 4
        return_metrics = false;
    end
    
    % Extract parameters
    K = x(1);
    a1 = x(2);
    b1 = x(3);
    a2 = x(4);
    b2 = x(5);
    a3 = x(6);
    b3 = x(7);
    zeta_n_ar = x(8);
    zeta_d_ar = x(9);
    zeta_n_r = x(10);
    zeta_d_r = x(11);
    
    % Define s for transfer functions
    s = tf('s');
    
    % Create the compensator
    % Lead compensators
    C1 = b1/a1 * (s + a1)/(s + b1);
    C2 = b2/a2 * (s + a2)/(s + b2);
    C3 = b3/a3 * (s + a3)/(s + b3);
    
    % Anti-notch filter for anti-resonance
    omega_ar = 4.601;
    notch_filter_ar = (s^2 + 2*zeta_n_ar*omega_ar*s + omega_ar^2) / (s^2 + 2*zeta_d_ar*omega_ar*s + omega_ar^2);
    
    % Notch filter for resonance
    omega_r = 8.347;
    notch_filter_r = (s^2 + 2*zeta_d_r*omega_r*s + omega_r^2) / (s^2 + 2*zeta_n_r*omega_r*s + omega_r^2);
    
    % Complete compensator
    C = K * C1 * C2 * C3 * notch_filter_ar * notch_filter_r;
    
    % Loop gain
    Lg = C * G;
    
    % Closed-loop system
    T = feedback(Lg, 1);
    
    % Compute bandwidth
    bw = bandwidth(T);
    
    % Compute margins
    [Gm, Pm, ~, ~] = margin(Lg);
    
    % Convert gain margin to dB
    gain_margin_dB = 20*log10(Gm);
    
    % Define target margins
    target_phase_margin = 40;  % degrees
    target_gain_margin = 10;   % dB

    % Calculate error metrics
    bandwidth_error = abs(bw - target_bandwidth);
    phase_margin_error = abs(Pm - target_phase_margin);
    gain_margin_error = abs(gain_margin_dB - target_gain_margin);

    % Total cost (weighted sum)
    cost = bandwidth_error + 100*phase_margin_error + 100*gain_margin_error;
    
    % Return metrics if requested
    if return_metrics
        metrics.bandwidth = bw;
        metrics.phase_margin = Pm;
        metrics.gain_margin = gain_margin_dB;
    end
end

%% Plotting Function
function plot_system(x, G, target_bandwidth)
    % Extract parameters
    K = x(1);
    a1 = x(2);
    b1 = x(3);
    a2 = x(4);
    b2 = x(5);
    a3 = x(6);
    b3 = x(7);
    zeta_n_ar = x(8);
    zeta_d_ar = x(9);
    zeta_n_r = x(10);
    zeta_d_r = x(11);
    
    % Define s for transfer functions
    s = tf('s');
    
    % Create the compensator
    % Lead compensators
    C1 = b1/a1 * (s + a1)/(s + b1);
    C2 = b2/a2 * (s + a2)/(s + b2);
    C3 = b3/a3 * (s + a3)/(s + b3);
    
    % Anti-notch filter for anti-resonance
    omega_ar = 4.601;
    notch_filter_ar = (s^2 + 2*zeta_n_ar*omega_ar*s + omega_ar^2) / (s^2 + 2*zeta_d_ar*omega_ar*s + omega_ar^2);
    
    % Notch filter for resonance
    omega_r = 8.347;
    notch_filter_r = (s^2 + 2*zeta_d_r*omega_r*s + omega_r^2) / (s^2 + 2*zeta_n_r*omega_r*s + omega_r^2);
    
    % Complete compensator
    C = K * C1 * C2 * C3 * notch_filter_ar * notch_filter_r;
    
    % Loop gain
    Lg = C * G;
    
    % Closed-loop system
    T = feedback(Lg, 1);
    
    % Calculate margins
    [Gm, Pm, Wcg, Wpc] = margin(Lg);
    gain_margin_dB = 20*log10(Gm);
    
    % Compute frequency responses
    w = logspace(-2, 2, 1000);
    [mag_lg, phase_lg] = bode(Lg, w);
    [mag_cl, phase_cl] = bode(T, w);
    mag_lg = squeeze(mag_lg);
    phase_lg = squeeze(phase_lg);
    mag_cl = squeeze(mag_cl);
    phase_cl = squeeze(phase_cl);
    
    % Calculate bandwidth
    bw = bandwidth(T);
    
    % Plot loop gain Bode
    figure('Name', 'Loop Gain Bode Plot');
    subplot(2,1,1);
    semilogx(w, 20*log10(mag_lg), 'LineWidth', 2);
    hold on;
    yline(0, 'r--', 'LineWidth', 1.5);
    if ~isempty(Wpc)
        semilogx(Wpc, 20*log10(mag_lg(find(w >= Wpc, 1))), 'ro', 'MarkerSize', 8);
        text(Wpc, 20*log10(mag_lg(find(w >= Wpc, 1))) + 5, ['GM = ' num2str(gain_margin_dB, '%.2f') ' dB']);
    end
    grid on;
    title('Loop Gain Magnitude');
    xlabel('Frequency (rad/s)');
    ylabel('Magnitude (dB)');
    
    subplot(2,1,2);
    semilogx(w, phase_lg, 'LineWidth', 2);
    hold on;
    yline(-180, 'r--', 'LineWidth', 1.5);
    if ~isempty(Wcg)
        idx = find(w >= Wcg, 1);
        if ~isempty(idx)
            semilogx(Wcg, phase_lg(idx), 'ro', 'MarkerSize', 8);
            text(Wcg, phase_lg(idx) + 10, ['PM = ' num2str(Pm, '%.2f') '°']);
        end
    end
    grid on;
    title('Loop Gain Phase');
    xlabel('Frequency (rad/s)');
    ylabel('Phase (deg)');
    
    % Plot closed-loop Bode
    figure('Name', 'Closed Loop Bode Plot');
    subplot(2,1,1);
    semilogx(w, 20*log10(mag_cl), 'LineWidth', 2);
    hold on;
    yline(-3, 'r--', 'LineWidth', 1.5);
    xline(bw, 'g--', 'LineWidth', 1.5);
    xline(target_bandwidth, 'b--', 'LineWidth', 1.5);
    grid on;
    title('Closed Loop Magnitude');
    xlabel('Frequency (rad/s)');
    ylabel('Magnitude (dB)');
    legend('Response', '-3 dB', ['Actual BW = ' num2str(bw/(2*pi), '%.2f') ' Hz'], ['Target BW = ' num2str(target_bandwidth/(2*pi), '%.2f') ' Hz']);
    
    subplot(2,1,2);
    semilogx(w, phase_cl, 'LineWidth', 2);
    grid on;
    title('Closed Loop Phase');
    xlabel('Frequency (rad/s)');
    ylabel('Phase (deg)');
    
    % Nyquist plot
    figure('Name', 'Nyquist Plot');
    nyquist(Lg);
    title('Nyquist Plot');
    axis equal;
    
    % Step response
    figure('Name', 'Step Response');
    step(T);
    grid on;
    title('Closed-Loop Step Response');
end 