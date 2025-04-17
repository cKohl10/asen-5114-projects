%% Sub Script for Problem 2 : Design State Variable Feedback Controller
% Authors: Thomas Dunnington, Owen Craig, Carson Kohbrenner
close all; clear; clc;

%% Extract State Space
SS = load("Data/state_space.mat");
A = SS.A';
B = SS.C';
C = SS.B';
D = SS.D';
ss_sys = ss(A,B,C,D);

% Initial guess: [Re, Im, p3,p4]
x0 = [-0.069, 4.6,-4,-3];

% Bounds: keep poles in stable LHP
lb = [-1.0, 3.5, -100, -100];
ub = [-0.01, 6.0, -0.05, -0.05];

plot_cost_space(x0, lb, ub, A, B, C);

function plot_cost_space(x0, lb, ub, A, B, C)
    % Create a grid of points in the x-space
    num_iter = 10;
    coords = {
        linspace(lb(1), ub(1), num_iter),
        linspace(lb(2), ub(2), num_iter),
        linspace(lb(3), ub(3), num_iter),
        linspace(lb(4), ub(4), num_iter),
    };
    % [X1, X2, X3, X4] = ndgrid(coords{1}, coords{2}, coords{3}, coords{4});
    % x_grid = [X1(:), X2(:), X3(:), X4(:)];

    % Calculate the cost for each point in the grid
    J_grid = zeros(num_iter, num_iter, num_iter, num_iter);
    max_count = num_iter^4;
    count = 0;
    for i = 1:num_iter
        for j = 1:num_iter
            for k = 1:num_iter
                for l = 1:num_iter
                    count = count+1;
                    if mod(count, round(max_count/10))==0
                        sprintf("%i / %i",count, max_count)
                    end
                    cost = pole_bandwidth_cost([coords{1}(i), coords{2}(j), coords{3}(k), coords{4}(l)], A, B, C);
                    if cost == 1e-6
                        cost = -1;
                    end
                    J_grid(i,j,k,l) = dB(cost);
                end
            end
        end
    end



    % Create a 4D matrix J from the grid and x_grid
    V = perms([1,2,3,4]);
    V = unique(V(:,1:2), 'rows');
    V = V(V(:,2) ~= 4,:);

    options = struct();
    options.dim_names = {'p1_{Re}', 'p1_{Im}', 'p3', 'p4'};
    options.const_idx = 1; % Use slice along Dim 4 (Time)
    options.frame_pause = 1;
    options.plot_type = 'imagesc'; % Use imagesc
    options.save_gif = true;  % Set to true to save as GIF
    % options.gif_filename = 'my_animation.gif';  % Optional: customize filename
    options.gif_delay = 1;  % Optional: set frame delay in seconds
    for i = 1:size(V,1)
        animate_slice_plot(J_grid, V(i,:), coords, options);
    end
        
end