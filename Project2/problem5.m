%% Main Script for Problem 5
% Authors: Thomas Dunnington, Owen Craig, Carson Kohbrenner
close all; clear; clc;


%% Simulate the full system with the observer
% Load the state space
load("Data/state_space.mat");

% Load the augmented state space
load('Data/observer_ss.mat')

% Load the controller
load('Data/closed_loop_ss.mat')

fig_size = [100, 100, 1000, 800];

% Step response
useSine = false;
u_amp = 0.1;
u_w = 0;
out = sim('Simulink/problem5_sim', 'StopTime', '10');
control_plot(out, fig_size, "Figures", "P5 Step Response");
observer_error_plot(out.t, out.x, out.xhat, "P5 Step Response");

% Step response 0.1 Hz
useSine = true;
u_w = 0.1 * 2 * pi; % rad/s
out = sim('Simulink/problem5_sim', 'StopTime', '50');
control_plot(out, fig_size, "Figures", "P5 Sine Response 0.1 Hz");
observer_error_plot(out.t, out.x, out.xhat, "P5 Sine Response 0.1 Hz");

% Step response 1 Hz
useSine = true;
u_amp = 60/10^(51/20);
u_w = 2 * pi; % rad/s
out = sim('Simulink/problem5_sim', 'StopTime', '20');
control_plot(out, fig_size, "Figures", "P5 Sine Response 1 Hz");
observer_error_plot(out.t, out.x, out.xhat, "P5 Sine Response 1 Hz");


