clear; close all; clc
rmse = @(org,est) norm(org(:) - est(:))./sqrt(numel(org));

%% create Test data 3D
rng('default');
[x,y,z] = meshgrid(linspace(-4,4,12));
A = sqrt(x.^2 + y.^2 + z.^2);
A(A<2) = 1; A(A>2) = 0;
sigma = 0.1;
W = sigma*randn(size(A));

Y = A + W;


figure(1), imagesc(A(:,:,6)), caxis([0 2]), colorbar
figure(2), imagesc(Y(:,:,6)), caxis([0 2]), colorbar

%% Display the 3D data

xs = -4:4;
ys = -4:4;
zs = -4:4;

figure(1), clf
contourslice(x,y,z,A, xs, ys, zs)
rotate3d on
view([-45 30])

figure(2), clf
contourslice(x,y,z,Y, xs(1:2:end), ys(1:2:end), zs(1:2:end))
rotate3d on
view([-45 30])
%% Denoise the 3d Data

blcksize = [12 12 2];
overlap = 6;
searchSize = [24 24 6];
threshold = 50;
lam = 0.01;
is2d = false;
iter = 10;
Est = X_small;
for i = 1:iter
    Est = Est + 0.1 * (X_small - Est);
    Est = lowRank3D(Est,blcksize,overlap,threshold,searchSize,is2d,lam);
    %fprintf('Denoised at step %d. RMSE = %1.4f \n', i, rmse(Est,A))
    i
end


%% 

figure(1), imagesc(Y(:,:,30)), caxis([0 2]), colorbar
figure(2), imagesc(Est(:,:,30)), caxis([0 2]), colorbar
figure(3), imagesc(A(:,:,30)), caxis([0 2]), colorbar

%% All slices at one voxel

figure(1), clf
subplot(2,2,1)
imagesc(reshape(A(30,:,:), [50 50]))
box off
colorbar
title('Clean Cube slice at x = 30')

subplot(2,2,2)
imagesc(reshape(A(:,30,:), [50 50]))
box off
colorbar
title('Clean Cube slice at y = 30')


subplot(2,2,3.5)
imagesc(reshape(A(:,:,30), [50 50]))
box off
colorbar
title('Clean Cube slice at z = 30')


figure(2), clf
subplot(2,2,1)
imagesc(reshape(Y(30,:,:), [50 50]))
box off
colorbar
title('Clean Cube slice at x = 30')

subplot(2,2,2)
imagesc(reshape(Y(:,30,:), [50 50]))
box off
colorbar
title('Clean Cube slice at y = 30')


subplot(2,2,3.5)
imagesc(reshape(Y(:,:,30), [50 50]))
box off
colorbar
title('Clean Cube slice at z = 30')


figure(3), clf
subplot(2,2,1)
imagesc(reshape(Est(30,:,:), [50 50]))
box off
colorbar
title('Clean Cube slice at x = 30')

subplot(2,2,2)
imagesc(reshape(Est(:,30,:), [50 50]))
box off
colorbar
title('Clean Cube slice at y = 30')


subplot(2,2,3.5)
imagesc(reshape(Est(:,:,30), [50 50]))
box off
colorbar
title('Clean Cube slice at z = 30')
