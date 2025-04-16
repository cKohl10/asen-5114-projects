%% Main Script for Problem 2 : Design State Variable Feedback Controller
% Authors: Thomas Dunnington, Owen Craig, Carson Kohbrenner
close all; clear; clc;

%% Extract State Space
SS = load("Data/state_space.mat");
A = SS.A';
B = SS.C';
C = SS.B';
D = SS.D';
ss_sys = ss(A,B,C,D);

%% Determine Desired Closed
% Desired closed-loop bandwidth = 1 Hz
j = sqrt(-1);
x_opt = optimize_bandwidth();
p1 = x_opt(1)+x_opt(2)*j;
p2 = x_opt(1)-x_opt(2)*j;
p3 = x_opt(3);     
p4 = x_opt(4);
Desired_poles = [p1, p2, p3, p4];

% Design Gain Matrix K to get desired poles
K = place(A,B,Desired_poles);
save('Data/K_matrix.mat', 'K');


%% Find Loop Gain, Closed Loop TF, and TF from Reference to Plant Input
% Loop Gain
s = tf('s');
%Lg_neg_test = K * inv(s * eye(size(A)) - A) * B;
Lg_neg = tf(ss(A, B, K, 0));

% Closed Loop TF
% F =  inv(C*inv(-A+B*K)*B);
% Cl_Tf = C*inv(s*eye(size(A))-A+B*K)*B*F;
Acl = A - B*K;
F = pinv(C * (Acl \ B));
Cl_Tf = tf(ss(Acl, B*F, C, 0));
        

% Reference to Plant Input
% Cl_u2r = inv(1+K*inv(s*eye(size(A))-A)*B)*F;
Cl_u2r = inv(1+Lg_neg)*F;
[mag_u2r, phase_u2r, wout_u2r] = bode(Cl_u2r);
mag_u2r = squeeze(mag_u2r);
phase_u2r = squeeze(phase_u2r);
input_limit_dB = 20*log10(67); 

%% Plot Bode Plot for Loop gain, Closed Loop, and Reference to Plant Input Response
%CL_bode_plot(Lg_neg, Cl_Tf ,'Figures/Problem2/', "Problem 2");
ss_bode_plots(A,B,C,D,K,F);

% Loop Gain Bode
% figure;
% bode(Lg_neg);
% grid on; % Turn on grid
% title('Bode Plot of Negative Loop Gain'); % Add a title
% xlabel('Frequency (rad/s)'); % X-axis label
% 
% % Closed Loop Bode
% figure;
% bode(Cl_Tf);
% grid on; % Turn on grid
% title('Bode Plot of Closed Loop TF'); % Add a title
% xlabel('Frequency (rad/s)'); % X-axis label

% Reference to Plant Input
% figure(6);
% set(gcf, 'Position', [100, 100, 500, 400]); 
% % Magnitude subplot
% subplot(2,1,1)
% semilogx(wout_u2r, db(mag_u2r), 'b', 'LineWidth', 2);
% hold on;
% xline(2*pi, 'k--', 'LineWidth', 1, 'Label', '1 Hz', 'LabelVerticalAlignment', 'bottom');
% yline(input_limit_dB, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Input Limit');
% ylabel('Magnitude (dB)');
% xlabel('Frequency (rad/s)');
% title('Reference to Plant Input Magnitude');
% legend('Cl\_u2r', 'Input Limit');
% grid on;
% xlim([wout_u2r(1), wout_u2r(end)]);
% % Phase subplot
% subplot(2,1,2)
% semilogx(wout_u2r, phase_u2r, 'b', 'LineWidth', 2);
% hold on;
% xline(2*pi, 'k--', 'LineWidth', 1, 'Label', '1 Hz', 'LabelVerticalAlignment', 'bottom');
% xlabel('Frequency (rad/s)');
% ylabel('Phase (deg)');
% title('Reference to Plant Input Phase');
% grid on;
% xlim([wout_u2r(1), wout_u2r(end)]);
% sgtitle('Bode Plot of Reference to Plant Input Transfer Function');

%% Plot Ol and Cl pole locations
% Calculate open-loop and closed-loop poles
open_loop_poles = pole(ss(A, B, C, 0));
closed_loop_poles = eig(A - B*K);

% Plot poles on complex plane
figure;
hold on;
grid on;
axis equal;
plot(real(open_loop_poles), imag(open_loop_poles), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot(real(closed_loop_poles), imag(closed_loop_poles), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Open-Loop and Closed-Loop Poles');
legend('Open-Loop Poles', 'Closed-Loop Poles');
xline(0, '--k');

%% Calculate the margins and bandwidth 
bw = bandwidth(Cl_Tf,-3.05);
[GM,PM] = margin(Lg_neg);
GM = 20*log10(GM);

%% Display State Variable Feedback Results
%clc;
fprintf('Closed-Loop Poles:\n');
for i = 1:length(closed_loop_poles)
    fprintf('%.4f %+.4fj\n', real(closed_loop_poles(i)), imag(closed_loop_poles(i)));
end
disp(['Gain Margin: ', num2str(GM), 'dB']);
disp(['Phase Margin: ', num2str(PM), 'deg']);
disp(['Closed-loop bandwidth: ', num2str(bw), ' rad/s']);
save('Data/closed_loop_ss.mat','F','K');

%% Functions 
function x_opt = optimize_bandwidth()
    % Load system
    SS = load("Data/state_space.mat");
    A = SS.A; B = SS.B; C = SS.C; D = SS.D;

    % Initial guess: [Re, Im, p3,p4]
    x0 = [-0.069, 4.06,-4,-3];
    
    % Bounds: keep poles in stable LHP
    lb = [-1000, 0.5,-40,-40];   % Avoid poles too close to jω
    ub = [-0.01, 20,-0.05,-0.05];

    % Run fmincon
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    x_opt = fmincon(@(x) pole_bandwidth_cost(x, A, B, C), x0, [], [], [], [], lb, ub, [], options);

    % % Show final results
    % fprintf('Optimal Pole Parts:\n');
    % disp(x_opt);
end
