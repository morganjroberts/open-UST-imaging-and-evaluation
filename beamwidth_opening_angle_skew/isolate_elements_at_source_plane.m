% DESCRIPTION:
%     Validation experiment to see which isolation technique provides the
%     best agreement with the measured plane in far field
%
% APPROACH:
%     1. write a function that builds the mask using known coordinates
%     2. write a function that finds the best mask position
%     3. write a function that applies the mask
%     4. re-project to the required distance
%     5. post-process the validation scan!
%     6. compare with the validation scan

close all
clearvars

[~, data_dir]    = getRepoDataPath();
scan_data_folder = [data_dir, filesep, 'field_scans'];

% -------------------------------------------------------------------------
% Load the datasets and backpropagate to a series of planes
filename       = 'field_scan_probe_A_all_channels_postprocessed';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
% M structure is the meausred pressure
M = load(input_filename, 'source_p', 'pressure', 'time_axis', 'Ny', 'Nx', 'Nt', 'source_z', 'c_water', ...
                                    'x_pos', 'y_pos', 'z_pos', 'dx', 'dt');

plotMeasurementPlane(M.source_p, M.x_pos, M.y_pos)
%%
[~,Ilin] = max(abs(M.source_p), [], 'all');
[iy, ix, it] = ind2sub(size(M.source_p), Ilin);

ix = 98;
truncatePressureInTime(M.source_p, M.time_axis, M.time_axis(end)*1e6, x=ix, y=iy, Clip=0.1);
ylim(gca, [-600,600]);

% UST probes: channel one has the highest (+ve x=coordinate)
%%
x8 = 0.0014; % estimated x-position of ch8 [m]
[~,i8] = findClosest(M.x_pos, x8);

ssp = sum(M.source_p.^2, 3);
max_p = max(M.source_p, [], 3);

figure;
plot(max_p(:, i8));

close all
% Build mask
xl = 144;
xr = 152;

lat_offset = -1;
xl = xl - lat_offset;
xr = xr + lat_offset;

y_offset = 1;
y_c = y_offset + round(M.Ny/2);
half_length = round(10e-3 / (2 * M.dx));
yt = y_c + half_length;
yb = y_c - half_length;

mask = zeros(M.Ny, M.Nx);
mask(yb:yt,xl:xr) = 1;
% mask(:,xl:xr) = 1;
mask_filt = imgaussfilt(mask, 0.6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preview
figure;
subplot(2, 2, 1);
imagesc(mask);
axis image
set(gca, 'YDir', 'normal')

subplot(2, 2, 3);
imagesc(mask_filt);
axis image
set(gca, 'YDir', 'normal')

subplot(2, 2, 2);
plot(mask(round(M.Ny/2),:))
hold on
plot(mask_filt(round(M.Ny/2),:))

subplot(2, 2, 4);
plot(mask(:, i8))
hold on
plot(mask_filt(:, i8))
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ssp slice before and after masking

figure;
plot(ssp(round(M.Ny/2),:));

% Preview the mask effect on SSP
mask_p = mask .* ssp;
mask_p_filt = mask_filt .* ssp;

% Preview
figure;
subplot(2, 2, 1);
imagesc(mask_p);
axis image
set(gca, 'YDir', 'normal')
title('Unfiltered mask');

subplot(2, 2, 3);
imagesc(mask_p_filt);
axis image
set(gca, 'YDir', 'normal')
title('Filtered mask');

subplot(2, 2, 2);
plot(mask_p(round(M.Ny/2),:))
hold on
plot(mask_p_filt(round(M.Ny/2),:))
plot(ssp(round(M.Ny/2),:))
xline(xl, 'k--');
xline(xr, 'k--');
legend({'Unfiltered mask', 'filtered mask', 'Unmasked'});

subplot(2, 2, 4);
plot(mask_p(:, i8,:))
hold on
plot(mask_p_filt(:, i8,:))
plot(ssp(:, i8,:))
xline(yb, 'k--');
xline(yt, 'k--');
legend({'Unfiltered mask', 'filtered mask', 'Unmasked'});

% Replicate mask for all time steps
mask_3d = repmat(mask_filt, [1, 1, M.Nt]);

% Apply mask to source plane
mask_source = mask_3d .* M.source_p;
%%
%%
% Project forward to validation plane

% load validation scan data
filename       = 'field_scan_probe_A_channel_8_validation_postprocessed';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
V = load(input_filename, 'pressure', 'z_pos', 'time_axis', 'Nt', 'Nx', 'Ny');

% NOTE: this does not need to be rounded to the nearest dx
z_proj = M.source_z + abs(V.z_pos) - abs(M.z_pos);

as_options = {'GridExpansion', 250, 'Plot', 0, 'Reverse', 0};
[~, val_p] = angularSpectrum(mask_source, M.dx, M.dt, z_proj, M.c_water, as_options{:});
%%
ssp_meas = sum(val_p.^2, 3);
ssp_val  = sum(V.pressure.^2, 3);

ssp_meas = ssp_meas / max(ssp_meas(:));
ssp_val = ssp_val / max(ssp_val(:));

figure;
subplot(3, 2, 1);
imagesc(ssp_meas);
axis image
colorbar
colormap(getBatlow);
title('Measured, reprojected');

subplot(3, 2, 3);
imagesc(ssp_val);
axis image
colorbar
colormap(getBatlow);
title('Validation');

diff = ssp_val - ssp_meas;
max_norm_err = 1e2* max((diff), [], 3) / max(ssp_val(:));

subplot(3, 2, 5);
imagesc(max_norm_err);
axis image
colorbar
colormap(getBatlow);
title('Difference [%]');

subplot(3, 2, 2);
hold on;
plot(ssp_val(round(M.Ny/2), :));
plot(ssp_meas(round(M.Ny/2), :));
legend({'Validation', 'Measured'});

subplot(3, 2, 4);
hold on;
plot(ssp_val(:, round(M.Nx/2)));
plot(ssp_meas(:, round(M.Nx/2)));
legend({'Validation', 'Measured'});


% 
% function buildElementMask(pitch, width, c0, options)
% % takes cartesian initial guess for c0
% % returns binary mask
% arguments
%     pitch
%     width
%     c0
%     options.Mode = 'X';
% end
% 
% switch options.Mode
%     case 'X'
% 
% 
% end