% FILENAME:
%     compare_reconstruction_and_phantom.m
% 
% DESCRIPTON:
%     Load a reconstruction result and compare it with the phantom.
%
% INPUT DATA FILENAMES:
%     <repo-dir>\image_reconstruction\sart_medfilt_11-Nov-2022-22-32-05_Nit-54_d-130mm.mat  reconstruction result object (instance of SartExperiment class)
%     <repo-dir>\phantom_design\phantom_sound_speed.mat                                     true sound speed map of phantom
%     
% FIGURE OUTPUT DIRECTORY:
%     <repo-dir>\image_reconstruction\figures
%
% DEPENDENCIES:
%     ust-sart https://github.com/ucl-bug/ust-sart
%     k-Wave   https://github.com/ucl-bug/k-wave
%
% ABOUT:
%     author:      - Morgan Roberts
%     last update: - 17/11/22

% Load the reconstruction result
sart_filename = 'sart_medfilt_11-Nov-2022-22-32-05_Nit-54_d-130mm.mat';
load(sart_filename)
sart = obj;
clear obj

% Load the phantom sound speed map
load('..\phantom_design\phantom_true_sound_speed.mat', 'phantom');

% Common scale limits for sound speed
cRange = [1420, 1580];

% Extract final reconstruction result
recon = sart.estimates(:,:,end);

% Make recon the same size as the phantom for better resolution, and so
% grid ticks are in exactly the same place
recon = resize(recon, size(phantom), 'bilinear');

fig = figure;

% Plot the reconstruction
subplot(1, 2, 1);
imagesc(1e3*sart.grid_x, 1e3*sart.grid_x, recon, cRange)
c = colorbar;
xlabel('x-position [mm]');
ylabel('y-position [mm]');
ylim(1e3*sart.grid_x([1, end]));
xlim(1e3*sart.grid_x([1, end]));
set(gca, 'YDir', 'normal');
colormap(getColorMap);
axis image
title('Reconstruction');

% Plot the phantom
subplot(1, 2, 2);
imagesc(1e3*sart.grid_x, 1e3*sart.grid_x, phantom, cRange)
c = colorbar;
xlabel('x-position [mm]');
ylabel('y-position [mm]');
ylim(1e3*sart.grid_x([1, end]));
xlim(1e3*sart.grid_x([1, end]));
ylabel(c, 'Sound Speed [m/s]');
set(gca, 'YDir', 'normal');
colormap(getColorMap);
axis image
title('True Values');

set(fig, 'Position', [404 278 932 602]);
set(fig,'renderer','Painters');
filename = [pwd, '\figures\phantom_reconstruction_comparison'];
print(fig, filename, '-dsvg');
print(fig, filename, '-depsc2');
