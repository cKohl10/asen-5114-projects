function [Lg, C, plt_fig] = loopgain(G,s, params, K, plot_flag)

    for i = 1:length(params)
        C_set(i) = comp_select(params{i}, s);
    end

    C = K;
    if plot_flag
        plt_fig = figure();
        title('Compensators Bode Plot');
    end

    for i = 1:length(C_set)
        C = C*C_set(i);
        
        if plot_flag
            [mag, phase, wout] = bode(C_set(i));
            phase = squeeze(phase);
            mag = squeeze(mag);
            subplot(2,1,1)
            semilogx(wout, db(mag), 'linewidth', 2, 'DisplayName', "Compensator " + num2str(i));
            hold on;
            title('Magnitude');
            xlabel('Frequency (rad/s)');
            ylabel('Amplitude (dB)');
            xlim([wout(1), wout(end)]);
            grid on;

            % Plot the phase plot
            subplot(2,1,2)
            semilogx(wout, phase, 'linewidth', 2, 'DisplayName', "Compensator " + num2str(i));
            hold on;
            title('Phase');
            xlabel('Frequency (rad/s)');
            ylabel('Phase (deg)');
            xlim([wout(1), wout(end)]);
            grid on;
        end
    end

    if plot_flag
        [mag, phase, wout] = bode(C);
        phase = squeeze(phase);
        mag = squeeze(mag);
        subplot(2,1,1)
        semilogx(wout, db(mag), 'k', 'linewidth', 2, 'DisplayName', 'Final Compensator');
        hold on;
        title('Magnitude');
        xlabel('Frequency (rad/s)');
        ylabel('Amplitude (dB)');
        xlim([wout(1), wout(end)]);
        grid on;
        legend('show');

        % Plot the phase plot
        subplot(2,1,2)
        semilogx(wout, phase, 'k', 'linewidth', 2, 'DisplayName', 'Final Compensator');
        hold on;
        title('Phase');
        xlabel('Frequency (rad/s)');
        ylabel('Phase (deg)');
        xlim([wout(1), wout(end)]);
        grid on;
        legend('show');
    end

    Lg = C*G;

end