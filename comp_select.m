function C = comp_select(params, s)
    % Select the compensator based on the type
    % type = 1: Lag Compensator / Lead Compensator
    % type = 2: Notch Filter
    % type = 3: Proportional Gain
    % type = 4: Integral Gain
    % type = 5: Derivative Gain


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