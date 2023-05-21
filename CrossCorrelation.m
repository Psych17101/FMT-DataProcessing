% Load the two images to be analyzed
im1 = imread('image1.tif');
im2 = imread('image2.tif');

% Define the size of the interrogation windows
win_size = 10; % in particles

% Loop over all interrogation windows
for i = 1:win_size:size(im1,1)-win_size
    for j = 1:win_size:size(im1,2)-win_size
        
        % Extract the current interrogation window from both images
        im1_win = im1(i:i+win_size-1, j:j+win_size-1);
        im2_win = im2(i:i+win_size-1, j:j+win_size-1);
        
        % Perform cross-correlation analysis to find the position of the correlation peak
        [x_offset, y_offset] = find_correlation(im1_win, im2_win); % write this function later
        
        % Use sub-pixel interpolation with a Gaussian fit to refine the peak position
        [x_offset, y_offset] = subpixel_interpolation(im1_win, im2_win, x_offset, y_offset); % write this function later
        
        % Calculate the velocity vector for this interrogation window
        [vx, vy] = calculate_velocity(x_offset, y_offset, win_size); % write this function later
        
        % Store the velocity vector for this interrogation window
        velocities(i:i+win_size-1, j:j+win_size-1, 1) = vx;
        velocities(i:i+win_size-1, j:j+win_size-1, 2) = vy;
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



function [x, y] =find_correlation(im1_win,im2_win)

end

function [x, y] = subpixel_interpolation(im1_win,im2_win,x_off,y_off)

end

function [x, y] = calculate_velocity(x_off,y_off,win_size)

end