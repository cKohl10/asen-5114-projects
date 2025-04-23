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

%% Load in Emperical Data
data = readmatrix('data/Spacecraft_spin_module_frequency_response_data.xlsx');
freq_exp = data(:,1)*2*pi;
amp = data(:,2);
emperical_mag = 20*log10(amp) - 20*log10(freq_exp);
emperical_phase = data(:,3) - pi/2;

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
T1 = tf(ss(A-L*C,B,K,0));
T2 = tf(ss(A-L*C,L,K,0));
P = tf(ss(A,B,C,0));
Lg_neg = (1/(1+T1))*T2*P;

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

%% Emperical Response
Temp = (1/(1+T1))*T2;
[mag_temp, phase_temp, wout] = bode(Temp,freq_exp);
phase_temp = squeeze(phase_temp);
mag_temp = squeeze(mag_temp);
emperical_lg_neg_mag = 10.^(emperical_mag./20).*mag_temp;
emperical_lg_neg_phase = emperical_phase + phase_temp;

%% Plotting
% Plot Bode
Lg_Cl_Bode_Plots(Lg_neg,Cl_Tf);
% Plot Nyquist
Nyquist_Plot(Lg_neg);
% Plot Step Response
StepResponsePlot(Cl_Tf, 'Step Response');
%Emperical Loop Gain Bode Plot 
emperical_bode_plot(emperical_lg_neg_mag,emperical_lg_neg_phase,wout);

