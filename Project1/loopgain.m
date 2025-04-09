function [Lg, C, plt_fig] = loopgain(G,s, params, K, plot_flag, fig_name)

    for i = 1:length(params)
        C_set(i) = comp_select(params{i}, s);
    end

    C = K;
    if plot_flag
        plt_fig = figure();
        title('Compensators Bode Plot');
    else
        plt_fig = [];
    end

    for i = 1:length(C_set)
        C = C*C_set(i);
        
        if plot_flag

            if params{i}(1) == 2 && params{i}(2) > params{i}(3)
                comp_name = "Anti-Notch Compensator";
            elseif params{i}(1) == 2 && params{i}(2) < params{i}(3)
                comp_name = "Notch Compensator";
            elseif params{i}(1) == 1 && params{i}(2) > params{i}(3)
                comp_name = "Lag Compensator";
            elseif params{i}(1) == 1 && params{i}(2) < params{i}(3)
                comp_name = "Lead Compensator";
            end

            [mag, phase, wout] = bode(C_set(i));
            phase = squeeze(phase);
            mag = squeeze(mag);
            subplot(2,1,1)
            semilogx(wout, db(mag), 'linewidth', 2, 'DisplayName', comp_name);
            hold on;
            title('Magnitude');
            xlabel('Frequency (rad/s)');
            ylabel('Amplitude (dB)');
            xlim([wout(1), wout(end)]);
            grid on;

            % Plot the phase plot
            subplot(2,1,2)
            semilogx(wout, phase, 'linewidth', 2, 'DisplayName', comp_name);
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
        legend('show', 'Location', 'eastoutside');

        % Plot the phase plot
        subplot(2,1,2)
        semilogx(wout, phase, 'k', 'linewidth', 2, 'DisplayName', 'Final Compensator');
        hold on;
        title('Phase');
        xlabel('Frequency (rad/s)');
        ylabel('Phase (deg)');
        xlim([wout(1), wout(end)]);
        grid on;
        legend('show', 'Location', 'eastoutside');
        saveas(plt_fig, fig_name);
    end

    Lg = C*G;

end