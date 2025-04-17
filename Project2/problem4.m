%% Main Script for Problem 4: Full state Observer
% Authors: Thomas Dunnington, Owen Craig, Carson Kohbrenner
close all; clear; clc;

%% Extract State Space
SS = load("Data/state_space.mat");
SVF = load("Data/closed_loop_ss.mat");
A = SS.A;
B = SS.B;
C = SS.C;
D = SS.D;
ss_sys = ss(A,B,C,D);
F = SVF.F;
K = SVF.K;

%% Calculate the Closed Loop and choose Observer Poles
Acl = A - B*K;
Cl_Tf = tf(ss(Acl, B*F, C, 0));
Cl_poles = pole(Cl_Tf);

obs_poles = 5* real(Cl_poles);
obs_poles(3) = obs_poles(3)+imag(Cl_poles(3))*1j;
obs_poles(4) = obs_poles(4)+imag(Cl_poles(4))*1j;
L = place(A', C', obs_poles)';

A_comb = [A-B*K , B*K; zeros(size(A)), A-L*C];
B_comb = [B*F;zeros(size(A,1),1)];
C_comb = [C,zeros(1,size(A,1))];
D_comb = D;

%% Analyse Step Response
ss_comb = ss(A_comb,B_comb,C_comb,D_comb);
Cl_Tf = tf(ss_comb);
ss_response = stepinfo(ss_comb,'SettlingTimeThreshold',0.05);


%% Calculate the Negative Loop Gain (Cut at Plant Input)
G1 = tf(ss(A-L*C,B,K,0));
G2 = tf(ss(A-L*C,L,K,0));
G3 = tf(ss(A,B,C,0));
Lg_neg = (1/(1+G1))*G2*G3;

%% Calculate the margins and bandwidth 
bw = bandwidth(Cl_Tf,-3.05);
[GM,PM] = margin(Lg_neg);
GM = 20*log10(GM);
closed_loop_poles = eig(A_comb);

%% Display Stability Margins and Bandwidth
if(ss_response.SettlingTime<3)
    disp(['The System meets the 5% Settling Time Requirement: ', num2str(ss_response.SettlingTime), ' sec'])
end

fprintf('Closed-Loop Poles:\n');
for i = 1:length(closed_loop_poles)
    fprintf('%.4f %+.4fj\n', real(closed_loop_poles(i)), imag(closed_loop_poles(i)));
end
disp(['Gain Margin: ', num2str(GM), 'dB']);
disp(['Phase Margin: ', num2str(PM), 'deg']);
disp(['Closed-loop bandwidth: ', num2str(bw), ' rad/s']);

%% Plotting
Lg_Cl_Bode_Plots(Lg_neg,Cl_Tf);
% figure;
% bode(Lg_neg);
% 
% figure;
% bode(Cl_Tf);

figure;
nyquist(Lg_neg);

figure; 
step(ss_comb);

