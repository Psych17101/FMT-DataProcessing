clear;
close all;
tic
im = imread('Pic1.tif');
im1 = im(1:end/2, :);
im1 = double(im1(:,:));
imwrite(im1, 'timeStep1.tif');
im2 = im(end/2+1:end, :);
im2 = double(im2(:,:));
imwrite(im2, 'timeStep2.tif');

im = imread('B00002_Calibration.tif');
im1 = im(1:end/2, :);
im1 = double(im1(:,:));
imwrite(im1, 'Calibration.tif');

