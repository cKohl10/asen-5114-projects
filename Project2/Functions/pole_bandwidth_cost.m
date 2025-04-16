function J = pole_bandwidth_cost(x, A, B, C)
    % Extract real and imaginary parts
    p1 = x(1) + 1j*x(2);
    p2 = x(1) - 1j*x(2);
    p3 = x(3);
    p4 = x(4);
    Desired_poles = [p1, p2, p3, p4];
    try
       % Place poles and compute closed-loop TF
        K = place(A, B, Desired_poles);
        Acl = A - B*K;
        F = pinv(C * (Acl \ B));
        Cl_Tf = tf(ss(Acl, B*F, C, 0));
        
        % Compute closed-loop bandwidth
        bw = bandwidth(Cl_Tf,-3);

        s = tf('s');
        % F = inv(C*inv(-A+B*K)*B);
        % Cl_Tf = C*inv(s*eye(size(A))-A+B*K)*B*F;
        Acl = A - B*K;
        sys_cl = ss(Acl, B, C, 0);
        F = 1 / dcgain(sys_cl);
        Cl_Tf = tf(sys_cl) * F;

        % Compute closed-loop bandwidth
        bw = bandwidth(Cl_Tf);
        % bw_penalty = abs(bw - 2*pi)^2;
        bw_penalty = max(0,bw-2*pi)^2;

        % Notch penalty
        [mag, ~] = bode(Cl_Tf, 4.6);
        mag_at_notch = squeeze(mag);  
        gain_dB = 20 * log10(mag_at_notch);
        notch_penalty = max(0, (-2 - gain_dB)^2);

        % w_notch = logspace(log10(3), log10(10), 300);  
        % [mag, ~] = bode(Cl_Tf, w_notch);
        % mag = squeeze(mag);
        % gain_dB = 20 * log10(mag);
        % drop_threshold = -2;
        % gain_drop = drop_threshold - gain_dB;
        % gain_drop(gain_drop < 0) = 0; 
        % notch_penalty = trapz(log10(w_notch), (gain_drop / 10).^2);


        % Compute gain at the notch frequency (4.6 rad/s)
        [mag, ~] = bode(Cl_Tf, 4.6);
        mag_at_notch = squeeze(mag);  
        gain_dB = 20 * log10(mag_at_notch);
        notch_penalty = max(0, (-2.5 - gain_dB)^2);

        % Penilize Small Margins
        sys_tf = ss(A, B, K, 0);
        Lg_neg = tf(sys_tf);  
        [GM,PM] = margin(Lg_neg);
        GM = 20*log10(GM);

        GM_target = 6;   
        PM_target = 45;

        gm_penalty = max(0, GM_target - GM)^2;
        pm_penalty = max(0, PM_target - PM)^2;

        % Create penality to attempt to keep the plant input at approx < 67 mNm
        Cl_u2r = inv(1+Lg_neg)*F;
        gain_limit = 67;
        w = logspace(-1, 2, 1000);
        [mag_u2r, ~] = bode(Cl_u2r, w);
        mag_u2r = squeeze(mag_u2r); 
        excess = mag_u2r - gain_limit;
        excess(excess < 0) = 0;
        input_penalty = trapz(w, log(1 + excess));

        % Final cost: bandwidth deviation + notch penalty
        J = 100*(bw - 2*pi)^2 + 200*notch_penalty + 0*input_penalty+ 2*gm_penalty + pm_penalty;
        
    catch
        J = 1e6;
    end
end




