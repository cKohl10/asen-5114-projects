function prob6(C, s)
    % Read the data from the file
    data = readmatrix('data/Spacecraft_spin_module_frequency_response_data.xlsx');
    freq = data(:,1)*2*pi;
    amp = data(:,2);
    amp_db_P = 20*log10(amp) - 20*log10(freq);
    amp_db_P0 = 20*log10(amp);
    phase_P = data(:,3) - pi/2;
    phase_P0 = data(:,3);

    % Transfer function with anti resonance, resonance pair
    DC_gain = 10^(-15/20);
    pole_1 = 0.3;
    pole_2 = 0.4;
    omega_ar1 = 4.601;   % anti-resonance frequency
    omega_r1 = 2;    % resonance frequency
    zeta_z = 0.015;
    zeta_p = 0.035;
    G_resonance = (s^2 + 2*zeta_z*omega_ar1*s + omega_ar1^2) / (s^2 + 2*zeta_p*omega_r1*s + omega_r1^2);
    G = DC_gain * G_resonance * 1 / (s + pole_1) * 1 / (s + pole_2);


    [mag, phase, wout] = bode(G, freq);
    phase = squeeze(phase);
    mag = squeeze(mag);

    % Determine Open loop Poles and Zeros
    OLPoles = pole(G);
    OLZeros = zero(G);

    % Calculate Gain and Phase Margin 
    Lg = -C*G;
    [GM,PM] = margin(-Lg); 
    GM = 20*log10(GM);
    Margin = [GM,PM]';

    % Determine Closed Loop Poles and Zeros
    clSys = feedback(-Lg, 1);  
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
    disp(' --- Problem 6 --- ')
    disp(['Gain Margin: ', num2str(GM), 'dB']);
    disp(['Phase Margin: ', num2str(PM), 'deg']);
    disp(['Closed-loop bandwidth: ', num2str(bw), ' rad/s']);

    % Create a Bode Plot for Loop Gain
    figure;
    hold on;
    bode(-Lg);
    title('P6 Loop Gain Bode Plot with Compensator');
    grid on;
    % saveas(gcf, 'figs/P3_BodePlot.png');

    % Create Bode plot of compesators and LG
    w = logspace(-1, 2, 1000);  % 0.1 to 100 rad/s
    figure;
    bode(C, w); 
    hold on;
    bode(-Lg,w)
    bode(G,w)
    legend('C','LG','Plant');
    title('P6 Bode Plots of Compensator, LG, Plant');
    grid on;


    % Create a Bode Plot for Loop Gain
    figure;
    hold on;
    bode(clSys);
    title('P6 Closed Loop Bode Plot with Compensator');
    grid on;
    % saveas(gcf, 'figs/P3_BodePlot.png');

    % Create a Nyquist Plot
    figure;
    hold on;
    nyquist(-Lg);
    title('P6 Loop Gain Nyquist Plot')
    grid on;

    % Plot Poles and Zeros 
    figure;
    hold on;
    plot(real(CLZeros), imag(CLZeros), 'ro', 'MarkerSize', 8);
    plot(real(CLPoles), imag(CLPoles), 'rx', 'MarkerSize', 10);
    plot(real(OLZeros), imag(OLZeros), 'bo', 'MarkerSize', 8);
    plot(real(OLPoles), imag(OLPoles), 'bx', 'MarkerSize', 10);

    xlabel('Real Axis');
    ylabel('Imaginary Axis');
    title('Pole-Zero Plot');
    grid on;
    axis equal; 
    legend('Closed Loop Zeros', 'Closed Loop Poles','Open Loop Zero','Open Loop Poles');
    hold off; 


    % Plot Tracking Error
    figure;
    semilogx(w, Ts_eval*100, 'LineWidth', 1.5);
    xlabel('Frequency (rad/s)');
    ylabel('Tracking Error (%)');
    title('Closed-Loop Tracking Error vs Frequency');
    grid on;
    % yline(5, '--r', '5% Error Threshold');
    ylim([0 100]);

end