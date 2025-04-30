%% Main Script for Problem 3
% Authors: Thomas Dunnington, Owen Craig, Carson Kohbrenner
close all; clear; clc;


%% Simulink Model %%
SS = load("Data/state_space.mat");
A = SS.A;
B = SS.B;
C = SS.C;
D = SS.D;
ss_sys = ss(A,B,C,D);

Cl = load("Data/closed_loop_ss.mat");
K = Cl.K;
F = Cl.F;

fig_size = [100, 100, 1000, 800];

% Step response
useSine = false;
u_amp = 0.1;
u_w = 0;
out = sim('model_response', 'StopTime', '10');
control_plot(out, fig_size, "Figures", "P3 Step Response");

% Step response 0.1 Hz
useSine = true;
u_w = 0.1 * 2 * pi; % rad/s
out = sim('model_response', 'StopTime', '50');
control_plot(out, fig_size, "Figures", "P3 Sine Response 0.1 Hz");

% Step response 1 Hz
useSine = true;
u_amp = 0.1;
u_w = 2 * pi; % rad/s
out = sim('model_response', 'StopTime', '50');
control_plot(out, fig_size, "Figures", "P3 Sine Response 1 Hz");


