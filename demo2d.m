%% 2D Image denoising using Enhanced Low-Rank Matrix Approximation
%
% Author: Ankit Parekh (ankit.parekh@nyu.edu)
% Last Edit: 11/4/16. 
%
% Please cite as: 
% A. Parekh and I. W. Selesnick. Enhanced Low-Rank Matrix Approximation. 
% IEEE Signal Processing Letters, 23(4):493-497, April 2016. 

clear; close all; clc
%% 2D example
rng('default');
S = double(imread('peppers.png'));
sigma = 30;
Y = S + sigma * randn(size(S));

figure(1)
subplot(1,2,1), imagesc(S), colormap(gray)
subplot(1,2,2), imagesc(Y), colormap(gray)
fprintf('Noisy Image PSNR = %2.2f dB \n', csnr(Y,S,0,0))
%%
blcksize = [6 6];
searchSize = [30 30];
overlap = 4;
threshold = 140;
lam = 140;
delta = 0.1;
is2d = true;
Est = Y;
for i = 1:10
    Est = Est + delta * (Y - Est);
    Est = lowRank3D(Est,blcksize,overlap,threshold,searchSize,is2d,lam);
    fprintf('Denoised at step %d. PSNR = %2.2f dB \n', i, csnr(Est,S,0,0))
end

%%
figure(1)
subplot(1,3,1), imshow(uint8(S))
subplot(1,3,2), imshow(uint8(Y))
subplot(1,3,3), imshow(uint8(Est))
