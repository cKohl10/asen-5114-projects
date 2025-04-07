function CL_bode_plot(Lg, save_location, prob_name)

    if ~isstruct(Lg)
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

    else
        mag_lg = Lg.lg_mag;
        phase_lg = Lg.lg_phase;
        wout_lg = Lg.lg_wout;
        mag_cl = Lg.cl_mag;
        phase_cl = Lg.cl_phase;
        wout_cl = Lg.cl_wout;

        % Calculate the closed loop bandwidth for the experimental data
        closed_loop_bandwidth = wout_cl(find(db(mag_cl) <= -3, 1, 'first')); %Bandwidth is the frequency where the magnitude is -3 dB
    end

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
    title(prob_name + ' Magnitude');
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
    title(prob_name + ' Phase');
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
    if save_location ~= ""
        saveas(gcf, save_location + prob_name + '_loop_gain_bode_plot.png');
    end


    %% Closed Loop Bode Plot
    figure;
    set(gcf, 'Position', [100, 100, 700, 500]); % Resize figure window
    subplot(2,1,1)
    semilogx(wout_cl, db(mag_cl), 'color', 'g', 'linewidth', 2);
    hold on;
    bandwidth_plot = xline(closed_loop_bandwidth, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--', 'label', closed_loop_bandwidth, 'LabelVerticalAlignment', 'bottom');
    yline(-3, 'linewidth', 1.5, 'color', 'r', 'linestyle', '--')
    title(prob_name + ' Magnitude');
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
    title(prob_name + ' Phase');
    xlabel('Frequency (rad/s)');
    ylabel('Phase (deg)');
    xlim([wout_cl(1), wout_cl(end)]);
    sgtitle('Closed Loop Bode Plot')
    grid on;
    if save_location ~= ""
        saveas(gcf, save_location + prob_name + '_closed_loop_bode_plot.png');
    end
end