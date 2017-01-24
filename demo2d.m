%% 2D Image denoising using Enhanced Low-Rank Matrix Approximation
%
% Author: Ankit Parekh (ankit.parekh@nyu.edu)
% Last Edit: 11/4/16. 
%
% Please cite as: 
% A. Parekh and I. W. Selesnick. Enhanced Low-Rank Matrix Approximation. 
% IEEE Signal Processing Letters, 23(4):493-497, April 2016. 

%% Initialize
clear; close all; clc
rng('default');
%% Load Test Image

S = double(imread('peppers.png'));
sigma = 30;
Y = S + sigma * randn(size(S));

figure(1)
subplot(1,2,1), imagesc(S), colormap(gray)
subplot(1,2,2), imagesc(Y), colormap(gray)
%% Run ELMA method for denoising Test Image

% Parameters 
blcksize = [8 8];
searchSize = [36 36];
overlap = 5;
threshold = 140;
lam = 160;
delta = 0.2;
is2d = true;
Est = Y;

% Main iterations
fprintf('Noisy Image PSNR = %2.2f dB \n', csnr(Y,S,0,0))
for i = 1:10
    Est = Est + delta * (Y - Est);
    Est = lowRank3D(Est,blcksize,overlap,threshold,searchSize,is2d,lam);
    fprintf('Denoised at step %d. PSNR = %2.2f dB \n', i, csnr(Est,S,0,0))
end

%% Plot the estimated image
figure(1)
subplot(1,3,1), imshow(uint8(S))
subplot(1,3,2), imshow(uint8(Y))
subplot(1,3,3), imshow(uint8(Est))
