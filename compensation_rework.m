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
data = readmatrix('data\Spacecraft_spin_module_frequency_response_data.xlsx');
freq_exp = data(:,1)*2*pi;
w_min = min(freq_exp);
% w_min = 0;
w_max = max(freq_exp);

figs = [figure, figure];
a_max = 0.1;
first_plot = true;
for a = 0.02:0.02:a_max
    for b = 1
        first_plot = false;
        transparency = a/a_max;
        params = {  % [1, a, 1];
                    % [1, 0.3, 20]; 
                    % [1, 0.4, 20]; 
                    [2, zeta_z, zeta_p, omega_r1, omega_ar1]
                }; 
        K = 100;
        Lg = loopgain(G,s, params, K);
        plot_lg(figs, Lg, w_min, w_max, transparency, first_plot);
    end
end

function plot_lg(figs, Lg, w_min, w_max, transparency, first_plot)
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
    if first_plot
        bandwidth_plot = xline(closed_loop_bandwidth, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--', 'label', closed_loop_bandwidth, 'LabelVerticalAlignment', 'bottom');
        yline(-3, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
        title('Magnitude');
        xlabel('Frequency (rad/s)');
        ylabel('Amplitude (dB)');
        xlim([wout_cl(1), wout_cl(end)]);
        legend([bandwidth_plot], 'Bandwidth');
        grid on;
    end

    % Plot the phase plot
    subplot(2,1,2)
    semilogx(wout_cl, phase_cl, 'color', [transparency, 0, 0], 'linewidth', 2);
    hold on;
    % xline(2*pi, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
    if first_plot
        title('Phase');
        xlabel('Frequency (rad/s)');
        ylabel('Phase (deg)');
        xlim([wout_cl(1), wout_cl(end)]);
        sgtitle('Closed Loop Bode Plot')
        grid on;
    end


    %% Nyquist Plot
    % Create a nyquist plot of the model
    % figure(figs(2));
    % nyquist(Lg, {w_min, w_max});
    % axis equal;


    % %% Notch Filter Plot
    % figure;
    % set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
    % subplot(2,1,1)
    % semilogx(wout_n, db(mag_n), 'b', 'linewidth', 2);
    % hold on;
    % title('Magnitude');
    % xlabel('Frequency (rad/s)');
    % ylabel('Amplitude (dB)');
    % xlim([wout_n(1), wout_n(end)]);
    % grid on;

    % Plot the phase plot
    % subplot(2,1,2)
    % semilogx(wout_n, phase_n, 'linewidth', 2);
    % hold on;
    % title('Phase');
    % xlabel('Frequency (rad/s)');
    % ylabel('Phase (deg)');
    % xlim([wout_n(1), wout_n(end)]);
    % sgtitle('Anti-Notch Filter Bode Plot')
    % grid on;
end

function Lg = loopgain(G,s, params, K)

    for i = 1:length(params)
        C_set(i) = comp_select(params{i}, s);
    end

    C = K;
    for i = 1:length(C_set)
        C = C*C_set(i);
    end

    % Negative Loop Gain
    Lg = C*G;

end

function C = comp_select(params, s)
    % Select the compensator based on the type
    % type = 1: Lag Compensator
    % type = 2: Lead Compensator
    % type = 3: Notch Filter

    C = 1;
    type = params(1);
    params = params(2:end);
    if type == 1
        % Lag Compensator (Pole is greater than zero)
        % a > b
        
        % Lead Compensator (Zero is greater than pole)
        % a < b
        a = params(1);
        b = params(2);
        C = (s/a + 1) / (s/b + 1);

    elseif type == 2
        % Notch Filter
        zeta_n = params(1);
        zeta_d = params(2);
        omega_n = params(3);
        omega_d = params(4);
        C = (s^2 + 2*zeta_n*omega_n*s + omega_n^2) / (s^2 + 2*zeta_d*omega_d*s + omega_d^2);

    elseif type == 3
        % Proportional Gain
        K = params(1);
        C = K;
    elseif type == 4
        % Integral Gain
        K = params(1);
        C = K / s;
    elseif type == 5
        % Derivative Gain
        K = params(1);
        C = K * s;
    end
end
