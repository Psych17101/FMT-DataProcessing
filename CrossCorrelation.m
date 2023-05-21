% Load the two images to be analyzed
im = imread('Pic1.tif');
im1 = im(1:end/2, :);
im2 = im(end/2+1:end, :);

figure(1)
image(im1)
figure(2)
image(im2)

% Define the size of the interrogation windows
win_size = 20; % in particles

iter = 1;
% Loop over all interrogation windows
for i = 1:win_size:size(im1,1)-win_size
    for j = 1:win_size:size(im1,2)-win_size
        
        % Extract the current interrogation window from both images
        im1_win = im1(i:i+win_size-1, j:j+win_size-1);
        im2_win = im2(i:i+win_size-1, j:j+win_size-1);

        correlation = xcorr2(im1_win, im2_win);
        l = 1:length(correlation);
        surf(l,l,correlation)
        
        % Perform cross-correlation analysis to find the position of the correlation peak
        [x_offset, y_offset] = find_correlation(im1_win, im2_win); % write this function later
        
        % Use sub-pixel interpolation with a Gaussian fit to refine the peak position
        [x_subpix, y_subpix] = subpixel_interpolation(im1_win, im2_win, x_offset, y_offset); % write this function later
        
        % Calculate the velocity vector for this interrogation window
        [u, v] = calculate_velocity(x_subpix, y_subpix, win_size); % write this function later
        
        % Store the velocity vector for this interrogation window
        velocities(i:i+win_size-1, j:j+win_size-1, 1) = u;
        velocities(i:i+win_size-1, j:j+win_size-1, 2) = v;
        iter = iter+1
    end
end

% Calculate the average velocity in each cell
avg_velocities = mean(velocities, [1 2]);

% Calculate the time between the two images
delta_t = 73e-6; % in seconds

% Scale the velocities by the time and the pixel size to get the wind velocity in m/s
pixel_size = 4.4e-6; % in meters
M = 0.05;
wind_velocities = avg_velocities .* pixel_size ./ delta_t;


% Define the grid for the flow field plot
[X, Y] = meshgrid(1:size(im1, 2), 1:size(im1, 1));

% Plot the flow field using quiver
figure;
quiver(X, Y, velocities(:,:,1), velocities(:,:,2));

% Set axis labels and title
xlabel('X');
ylabel('Y');
title('Flow Field');

% Adjust the axis limits if needed
% xlim([xmin, xmax]);
% ylim([ymin, ymax]);

function [x, y] = find_correlation(im1_win,im2_win)
    % Perform cross-correlation between the two images
    correlation = xcorr2(im1_win, im2_win);
    % Find the position of the correlation peak
    [peak_value, peak_index] = max(correlation(:));
    [y_peak, x_peak] = ind2sub(size(correlation), peak_index);

    % Calculate the offset from the center of the interrogation window
    x = x_peak - size(im1_win,2);
    y = y_peak - size(im1_win,1);
end

function [x_offset, y_offset] = subpixel_interpolation(im1, im2, x_peak, y_peak)
    % Define the size of the search area around the correlation peak
    search_size = 3;

    % Determine the valid range for the subregion extraction
    [im_height, im_width] = size(im1);
    min_y = max(1, y_peak - search_size);
    max_y = min(im_height, y_peak + search_size);
    min_x = max(1, x_peak - search_size);
    max_x = min(im_width, x_peak + search_size);

    % Extract the subregion around the correlation peak
    subregion_im1 = im1(min_y:max_y, min_x:max_x);
    subregion_im2 = im2(min_y:max_y, min_x:max_x);

    % Calculate the center of mass for the correlation peak
    [peak_y, peak_x] = find(subregion_im2 == max(subregion_im2(:)));
    peak_y = peak_y(1);
    peak_x = peak_x(1);

    % Calculate the subpixel offset using weighted average
    weight_sum = sum(subregion_im2(:));
    x_offset = peak_x - (search_size + 1) + sum(double(1:size(subregion_im2, 2)) .* double(subregion_im2(peak_y, :))) / weight_sum;
    y_offset = peak_y - (search_size + 1) + sum(double(1:size(subregion_im2, 1))' .* double(subregion_im2(:, peak_x))) / weight_sum;

    % Adjust the offsets to be relative to the full image coordinates
    x_offset = x_peak + x_offset - 1;
    y_offset = y_peak + y_offset - 1;
end

function [u, v] = calculate_velocity(x_subpix,y_subpix,win_size)
pitch = 4e-6;
M = 0.05;
t = 74e-3; % may not be right
u = x_subpix.*pitch./(M.*t);
v = y_subpix.*pitch./(M.*t);

end