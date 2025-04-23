%% Main Script for Problem 2 : Design State Variable Feedback Controller
% Authors: Thomas Dunnington, Owen Craig, Carson Kohbrenner
close all; clear; clc;

%% Extract State Space
load("Data/state_space.mat");
ss_sys = ss(A,B,C,D);

%% Determine Closed and Open Loop Poles 
OlPoles = pole(ss_sys);
wn = 4.601;
zeta = 0.015;
Desired_poles = [-zeta*wn+wn*sqrt(1-zeta^2)*j,-zeta*wn-wn*sqrt(1-zeta^2)*j,-4,-3];


%% Design Gain Matrix K to get desired poles
K = place(A,B,Desired_poles);
F = inv(C * inv(-A + B*K) * B);


%% Plot the Open Loop and Closed Loop Frequency Responses
% Plot the loop gain and closed loop
ss_bode_plots(A, B, C, D, K, F);


%% Find Loop Gain and Closed Loop TF
s = tf('s');
Cl = C * inv(eye(size(A)).*s - A + B*K) * B * F;

% Calculate the bandwidth of the closed loop transfer function
Cl_bandwidth = bandwidth(Cl)