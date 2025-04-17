function animate_slice_plot(J, V, dim_ticks, options)
    J_r = zeros(length(J), length(J), length(J));
    J_h = zeros(length(J), length(J));
    
    % Determine coordinates based on first dimension we're slicing through
    switch V(1)
        case 4
            J_r(:,:,:) = J(:,:,:,round(end/2));
            coords = options.dim_names([1,2,3]);
            dim_ticks = dim_ticks([1,2,3]);
        case 3
            J_r(:,:,:) = J(:,:,round(end/2),:);
            coords = options.dim_names([1,2,4]);
            dim_ticks = dim_ticks([1,2,4]);
        case 2
            J_r(:,:,:) = J(:,round(end/2),:,:);
            coords = options.dim_names([1,3,4]);
            dim_ticks = dim_ticks([1,3,4]);
        case 1
            J_r(:,:,:) = J(round(end/2),:,:,:);
            coords = options.dim_names([2,3,4]);
            dim_ticks = dim_ticks([2,3,4]);
    end

    f = figure('Position', [100, 100, 800, 600]);
    h = heatmap(ones(size(J_h)), 'ColorbarVisible', 'on');
    colormap(jet);
    % Set the title through the heatmap object instead of using the title function
    h.Title = 'Cost Animation';
    h.FontSize = 14;
    
    % Setup for gif saving
    if isfield(options, 'save_gif') && options.save_gif
        gif_filename = sprintf('Figures/J_%i_over_%s.gif', length(dim_ticks), coords{V(2)});
        if isfield(options, 'gif_filename')
            gif_filename = options.gif_filename;
        end
        frame_delay = 0.1; % Default delay between frames in seconds
        if isfield(options, 'gif_delay')
            frame_delay = options.gif_delay;
        end
    end
    
    % Animation loop
    for i = 1:size(J_r, 1)
        % Extract 2D slice for current frame
        switch V(2)
            case 3
                J_h(:,:) = J_r(:,:,i);
                local_coords = {coords{1}, coords{2}};
                local_ticks = dim_ticks([1,2]);
                title_str = sprintf('%s = %f', coords{3}, dim_ticks{3}(i));
            case 2
                J_h(:,:) = J_r(:,i,:);
                local_coords = {coords{1}, coords{3}};
                local_ticks = dim_ticks([1,3]);
                title_str = sprintf('%s = %f', coords{2}, dim_ticks{2}(i));
            case 1
                J_h(:,:) = J_r(i,:,:);
                local_coords = {coords{2}, coords{3}};
                local_ticks = dim_ticks([2,3]);
                title_str = sprintf('%s = %f', coords{1}, dim_ticks{1}(i));
        end
        
        % Update heatmap with current data
        h.ColorData = J_h;
        h.XLabel = local_coords{1};
        h.YLabel = local_coords{2};
        h.XData = local_ticks{1};
        h.YData = local_ticks{2};
        h.Title = sprintf('%s (Slice %d of %d)', title_str, i, size(J_r, 1));
        globalMin = min(J_r(:));
        globalMax = max(J_r(:));
        h.ColorLimits = [globalMin, globalMax];
        
        drawnow;
        
        % Capture frame for GIF
        if isfield(options, 'save_gif') && options.save_gif
            frame = getframe(f);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);
            
            % Write to GIF file
            if i == 1
                imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', frame_delay);
            else
                imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
            end
        end
        
        pause(options.frame_pause);
    end
end