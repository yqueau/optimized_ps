%Clear workspace
clear
close all
%~ clc

%% Settings

dataset = 'Beethoven';

% Classic PS settings
lambda_classical = 1e-6; % Prior weight for integration

% Nonconvex postprocessing settings
lambda_nonconvex = 1e-6; % Prior weight, smaller values result in less influence of 
                         % the prior from classic PS
tol_global = 1e-8;       % Tolerance for the outer, alternating optimisation
tol_iPiano = 1e-8;       % Tolerance for the inner optimisation of the iPiano method
maxit_alternate = 1000;    % Max. number of global iterations

%% Computation
% Load PS functions
addpath Toolbox/

% Load data :
% I : nrows x ncols x nimgs   === Tensor of images - I(:,:,i) is image #i
% S : nimgs x 3               === Light matrix - S(i,:) = [si_x,si_y,zi_z] is light #i
%                                 Axis O--->y  (ij-coordinates)
%                                      |
%                                      x
% mask : nrows x ncols        === Binary mask of the object
load(['Datasets/',dataset,'.mat']);
I = I/255;
imask = find(mask>0);

% Make finite difference coefficient matrix
% M : 2.npix x npix (sparse)  === such that Dx(z) = M(1:2:end-1,:)*z(find(mask))
%                                       and Dy(z) = M(2:2:end,:)*z(find(mask))
M = make_gradient(mask);

% Classic PS
z0 = zeros(size(mask)); % Prior for integration
[z_classical,rho_classical,N_nonintegrable] = classical_ps(I,S,mask,lambda_classical,z0,M);

% Refined PS
[z_nonconvex,rho_nonconvex] = refine_PS(I,S,mask,lambda_nonconvex,z_classical,rho_classical,M,tol_global,tol_iPiano,maxit_alternate);

%% Visualisation

% Normals from depth map, classic PS model
Dz = M*z_classical(imask);
zx = zeros(size(mask));
zy = zeros(size(mask));
zx(imask) = Dz(1:2:end-1);
zy(imask) = Dz(2:2:end);

% Reprojected error
denom = sqrt(1+zx.^2+zy.^2);
N_classical = cat(3,-zx./denom,-zy./denom,1./denom);
SE_classical =  reprojection_PS(S,rho_classical,N_classical,I,imask);
MSE_classical = mean(SE_classical(imask));
SE_classical(mask==0) = Inf;
Reproj_classical = sum(SE_classical(imask));

figure;
imagesc(SE_classical,[0 0.0025])
colorbar
colormap gray
title(sprintf('Classic PS model, MSE = %.3e',MSE_classical))
drawnow

% Normals from depth map, nonconvex PS model
Dz = M*z_nonconvex(imask);
zx = zeros(size(mask));
zy = zeros(size(mask));
zx(imask) = Dz(1:2:end-1);
zy(imask) = Dz(2:2:end);

% Reprojected error
denom = sqrt(1+zx.^2+zy.^2);
N_nonconvex = cat(3,-zx./denom,-zy./denom,1./denom);
SE_nonconvex =  reprojection_PS(S,rho_nonconvex,N_nonconvex,I,imask);
MSE_nonconvex = mean(SE_nonconvex(imask));
SE_nonconvex(mask==0) = Inf;
Reproj_nonconvex = sum(SE_nonconvex(imask));

figure;
imagesc(SE_nonconvex,[0 0.0025])
colorbar
colormap gray
title(sprintf('Nonconvex PS model, MSE = %.3e',MSE_nonconvex))

% Shape, nonconvex PS model
figure
surfl(z_classical,[0 90]);axis equal;axis off; axis ij; shading flat; colormap gray;view(-35,40);
title('Classic PS model, shape');

% Shape, nonconvex PS model
figure
surfl(z_nonconvex,[0 90]);axis equal;axis off; axis ij; shading flat; colormap gray;view(-35,40);
title('Nonconvex PS model, shape');
drawnow
