%% Main Script for Problem 2 : Design State Variable Feedback Controller
% Authors: Thomas Dunnington, Owen Craig, Carson Kohbrenner
close all; clear; clc;

%% Extract State Space
SS = load("Data/state_space.mat");
A = SS.A;
B = SS.B;
C = SS.C;
D = SS.D;
ss_sys = ss(A,B,C,D);

%% Determine Closed and Open Loop Poles 
j = sqrt(-1);
OlPoles = pole(ss_sys);
wn = 2*pi;
zeta = 0.1;
Desired_poles = [-zeta*wn+wn*sqrt(1-zeta^2)*j,-zeta*wn-wn*sqrt(1-zeta^2)*j,-4,-3];

%% Design Gain Matrix K to get desired poles
K = place(A,B,Desired_poles);
F = 1;

%% Find Loop Gain and Closed Loop TF
s = tf('s');
sys_tf = ss(A, B, K, 0);
Lg_neg = tf(sys_tf);  

Acl = A - B*K;
Bcl = B*F;
Cl_Tf = tf(ss(Acl, Bcl, C, 0));

%% Plot Bode Plot for Loop gain and Closed Loop
figure;
bode(Lg_neg);
grid on; % Turn on grid
title('Bode Plot of Lg\_neg'); % Add a title
xlabel('Frequency (rad/s)'); % X-axis label

figure;
bode(Cl_Tf);
grid on; % Turn on grid
title('Bode Plot of Cl\_Tf'); % Add a title
xlabel('Frequency (rad/s)'); % X-axis label

Cl_bandwidth = bandwidth(Cl_Tf)