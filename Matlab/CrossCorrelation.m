% Load the two images to be analyzed
clear;
close all;
tic
im = imread('Pic1.tif');
im1 = im(1:end/2, :);
im1 = double(im1(:,:));
%im1 = flip(a1',1);
im2 = im(end/2+1:end, :);
im2 = double(im2(:,:));
%im2 = flip(a2',1);

figure(1)
image(im1)
figure(2)
image(im2)

out_name = 'pic1_converted.png';   %the name of the desired output
imwrite(im1, out_name);  %write it out
out_name = 'pic2_converted.png';   %the name of the desired output
imwrite(im2, out_name);  %write it out

% Define the size of the interrogation windows
win_size = 32; % pixels


scale = 0.5;
tscale = 73e-6;

itery = 1;
iterx = 1;
perc = 0;
fprintf('Entering Cycle \n')
% Loop over all interrogation windows
for i = 1:win_size/2:size(im1,1)-win_size % y axis
    iterx ;
    for j = 1:win_size/2:size(im1,2)-win_size %x axis 
        itery;
        % Extract the current interrogation window from both images
        im1_win = im1(i:i+win_size-1, j:j+win_size-1);
        im2_win = im2(i:i+win_size-1, j:j+win_size-1);
        int_center_x = (i+i+win_size-1)/2;
        int_center_y = (j+j+win_size-1)/2;
        a_win = im1_win-mean(im1_win(:));
        b_win = im2_win-mean(im2_win(:));
        % Perform cross-correlation analysis to find the position of the correlation peak
        [x_offset, y_offset,x_peak,y_peak,correlation,peak_index] = find_correlation(a_win, b_win); % write this function later

        % Use sub-pixel interpolation with a Gaussian fit to refine the peak position
        [x_subpix,y_subpix] = subpixel_int(correlation,x_peak,y_peak,x_offset,y_offset); % write this function later
        
        % Calculate the velocity vector for this interrogation window
        vecx(i:i+win_size-1, j:j+win_size-1, 1) = x_subpix;
        vecy(i:i+win_size-1, j:j+win_size-1, 1) = y_subpix;
        vec(i:i+win_size-1, j:j+win_size-1, 1) = x_subpix^2+y_subpix^2;
        
        [u, v] = calculate_velocities(x_subpix, y_subpix, win_size); % write this function later
      
        % Store the velocity vector for this interrogation window
        velocities(i:i+win_size-1, j:j+win_size-1, 1) = u;
        velocities(i:i+win_size-1, j:j+win_size-1, 2) = v;
        vel(i:i+win_size-1, j:j+win_size-1) =vec(i:i+win_size-1, j:j+win_size-1)*scale/tscale;

        itery = itery + 1;
        p = (length(im1_win)/2)*(length(im1_win)/2)/(length(im1(:,1))*length(im1(1,:)));
        perc = perc + p;
    end

    disp([num2str(floor(100*(perc))),'% passed']); %Progress (%)
    iterx = iterx + 1;
    itery = 1;
end
fprintf('Done!')
%%

%%
% Define the grid for the flow field plot
[X, Y] = meshgrid(1:length(velocities(:,:,1)), 1:length(velocities(:,1,:)));

% Plot the flow field using quiver
figure;
quiver(X, Y, velocities(:,:,1), velocities(:,:,2));

% Set axis labels and title
xlabel('X');
ylabel('Y');
title('Flow Field');


% Adjust the axis limits if needed
% xlim([0, i_im]);
% ylim([0, j]);

% Calculate the magnitude of velocities
magnitude = sqrt(velocities(:,:,1).^2 + velocities(:,:,2).^2);

threshold = 15; 
high_values = magnitude > threshold;
magnitude(high_values) = threshold;

[m, n] = size(magnitude);
magnitude_flipped = flip(magnitude(1:m, 1:n));

% Create a grid for the flow field plot
[X, Y] = meshgrid(1:size(velocities,2), 1:size(velocities,1));

% Create a figure and plot the heat map
figure;
contourf(X, Y, magnitude_flipped);
colormap(jet);
colorbar; % Add a colorbar to indicate velocity magnitude


% Set axis labels and title
xlabel('X');
ylabel('Y');
title('Velocity Magnitude Heat Map');
toc

%% Functions
function [x, y,x_peak,y_peak,correlation,peak_index] = find_correlation(im1_win,im2_win)
    % Perform cross-correlation between the two images
    correlation = normxcorr2(im1_win, im2_win);
    % Find the position of the correlation peak
    [peak_value, peak_index] = max(correlation(:));
    [x_peak, y_peak] = ind2sub(size(correlation), peak_index);

    % Calculate the offset from the center of the interrogation window
    x = x_peak - size(im1_win,1);
    y = y_peak - size(im1_win,2);
end

function [x_subpix, y_subpix] = subpixel_int(correlation, x_peak, y_peak,x_offset,y_offset)
maxi = correlation(x_peak, y_peak);   
% Conditions on x
        if (x_peak < length(correlation) && x_peak < length(correlation)-1)&& x_peak > 2 
            % Analysis if x_peak is within bounds of cross correlation.
            maxi_0x = correlation(x_peak - 2, y_peak);
            maxi_1x = correlation(x_peak - 1, y_peak);
            maxi_2x = correlation(x_peak + 1, y_peak);
            maxi_3x = correlation(x_peak + 2, y_peak);
        elseif x_peak == length(correlation)
            % if peak is at the right edge of cross-correlation window
            maxi_0x = correlation(x_peak - 2, y_peak);
            maxi_1x = correlation(x_peak - 1, y_peak);
            maxi_2x = 0;
            maxi_3x = 0;
        elseif x_peak < 2 && x_peak ~= 1 % Might be the wrong ones
            % if peak is on the right end of cross-correlation window 
            maxi_0x = 0;
            maxi_1x = correlation(x_peak - 1, y_peak);
            maxi_2x = 0;
            maxi_3x = 0;
        else % Should never be the case - but put just in case
            maxi_0x = 0;
            maxi_1x = 0;
            maxi_2x = 0;
            maxi_3x = 0;
        end
% Conditions on y
        if (y_peak < length(correlation) && y_peak < length(correlation)-1) && y_peak > 2
            % Analysis if x_peak is within bounds of cross correlation.
            maxi_0y = correlation(x_peak, y_peak - 2);
            maxi_1y = correlation(x_peak, y_peak - 1);
            maxi_2y = correlation(x_peak, y_peak + 1);
            maxi_3y = correlation(x_peak, y_peak + 2);
        elseif y_peak == length(correlation)
            % if peak is at the right edge of cross-correlation window
            maxi_0y = correlation(x_peak, y_peak - 2);
            maxi_1y = correlation(x_peak, y_peak - 1);
            maxi_2y = 0;
            maxi_3y = 0;
         elseif y_peak < 2 && y_peak ~= 1
            % if peak is on the right end of cross-correlation window 
            maxi_0y = 0;
            maxi_1y = correlation(x_peak, y_peak - 1);
            maxi_2y = 0;
            maxi_3y = 0;
        else % Should never be the case - but put just in case
            maxi_0y = 0;
            maxi_1y = 0;
            maxi_2y = 0;
            maxi_3y = 0;
        end

    x = 0:4;
    z_1 = double([maxi_0x , maxi_1x, maxi, maxi_2x, maxi_3x]);
    y = [y_peak -2 ,y_peak - 1, y_peak, y_peak + 1,y_peak + 2];
    y = 0:4;
    z_2 = [maxi_0y,maxi_1y, maxi, maxi_2y,maxi_3y];

    gaussEqn = 'a*exp(-((x-b)/c)^2)';
    fitresult = fit(x',z_1',gaussEqn,'StartPoint', [0.00001, 1, 1],'Upper',[maxi,10,10]);
    % Extract the fitted parameters
    A_x = fitresult.a;
    mu_x = fitresult.b;
    sigma_x = fitresult.c;
    x_subpix = x_offset + mu_x - 1;


    % Plot the data and the fitted curve
%     figure(7)
%     plot(x, z_1, 'o', 'DisplayName', 'Data');
%     hold on;
%     x_fit = linspace(0, 4, 1000);
%     y_fit = feval(fitresult, x_fit);
%     plot(x_fit, y_fit, 'r', 'DisplayName', 'Gaussian Fit');
%     legend('Location', 'best');
%     xlabel('x');
%     ylabel('z_1');
%     title('Gaussian Curve Fitting');
%     hold off;

   
    fitresult = fit(y',z_2',gaussEqn,'StartPoint', [0.000001, 1, 1],'Upper',[maxi,2,2]);
    % Extract the fitted parameters
    A_y = fitresult.a;
    mu_y = fitresult.b;
    sigma_y = fitresult.c;
    y_subpix = y_offset + mu_y - 1;

%  % Plot the data and the fitted curve
    
%    figure(8)
%    plot(y, z_2, 'o', 'DisplayName', 'Data');
%    hold on;
%    x_fit = linspace(0, 4, 1000);
%    y_fit = feval(fitresult, x_fit);
%    plot(x_fit, y_fit, 'r', 'DisplayName', 'Gaussian Fit');
%    legend('Location', 'best');
%    xlabel('x');
%    ylabel('z_1');
%    title('Gaussian Curve Fitting');
%    hold off;

   

end


function [u, v] = calculate_velocities(x_subpix,y_subpix,win_size)
pitch = 4.4e-6;
M = 0.05;
t = 73e-6; % may not be right
u = x_subpix.*pitch./(M.*t);
v = y_subpix.*pitch./(M.*t);

end