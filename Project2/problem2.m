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
zeta = 0.7;
Desired_poles = [-zeta*wn+wn*sqrt(1-zeta^2)*j,-0.4-6*j,-2.4,-2.3];

%% Design Gain Matrix K to get desired poles A-BK+BFR
K = place(A,B,Desired_poles);
s = tf('s');
Lg_neg  = K*inv((s*eye(size(A))-A))*B;

%Cl_Tf = C*inv(s*eye(size(A))-A+B*K)*B*F*R;

figure;
bode(Lg_neg);


