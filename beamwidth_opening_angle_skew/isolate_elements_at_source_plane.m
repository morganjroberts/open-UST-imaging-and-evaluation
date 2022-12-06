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

% UST probes: channel one has the highest (+ve x=coordinate)
%%

% Estimate where element 8 is (coordinates extracted visually from ssp)
x0 = 0.0014;  
y0 = 0.00035;
[~,ix0] = findClosest(M.x_pos, x0);
[~,iy0] = findClosest(M.y_pos, y0);
% ix0 = ix0+1;
% iy0 = iy0-2;
% Build a mask to isolate element 8

% Create ssp matrix for fine mask alignment
ssp = sum(M.source_p.^2, 3);

% remove border
x_border = 80;
y_border = 10;
ssp_cut = ssp( y_border+1:end-y_border, x_border+1:end-x_border );

% upsample
upsample_factor = 7; % MUST BE ODD
sz              = size(ssp_cut);
sz_up           = sz * upsample_factor;
ssp_up          = interpftn(ssp_cut, sz_up);

% figure;
% subplot(3, 1, 1);
% imagesc(ssp_cut);
% axis image
% 
% subplot(3, 1, 2);
% imagesc(ssp_up);
% axis image
% 
% subplot(3, 1, 3);
% hold on;
% plot(1:upsample_factor:upsample_factor*size(ssp_cut, 2), ssp_cut(iy0,:), 'kx');
% plot(ssp_up( ( (iy0 - 1) * upsample_factor) + 1,: ), 'b')
%-------------------------------------------------------------------------
% Visual identification of element centre (ch8) after border cutting
ix0 = 68; %68
iy0 = 21; %21


ix0_up = ix0 * upsample_factor;
iy0_up = iy0 * upsample_factor;

dx_up = M.dx / upsample_factor;
Nx_up = size(ssp_cut, 2) * upsample_factor;
Ny_up = size(ssp_cut, 1) * upsample_factor;

% source dimesions must be multiple of dx_up
mask_up = createElementMask(Nx_up, Ny_up, dx_up, 0.9e-3, 9e-3, ix0_up, iy0_up);


mask = imresize(mask_up, sz, 'bilinear', 'Antialiasing', true);
mask_f = interpftn(mask_up, sz);

figure;
subplot(3, 1, 1);
imagesc(mask_up);
axis image

subplot(3, 1, 2);
imagesc(mask);
axis image

i_vec = (1:upsample_factor:Nx_up) + ((upsample_factor - 1) / 2);

subplot(3, 1, 3);
hold on;
plot(i_vec, mask(iy0,:), 'k');
plot(mask_up(iy0_up,: ), 'b');


max_shift = 15;
shifts = -max_shift:max_shift;
Npos = length(shifts);
r = zeros(Npos, Npos);
for xdx = 1:Npos
    for ydx = 1:Npos
        x_shift = shifts(xdx);
        y_shift = shifts(ydx);
        mask_shift = circshift(mask_up, [y_shift, x_shift]);
        r(ydx, xdx) = sum(mask_shift .* ssp_up, 'all');
    end
end

[~,I] = max(r, [], 'all');
[y_opt, x_opt] = ind2sub(size(r), I);
ix = shifts(x_opt);
iy = shifts(y_opt);

figure;
imagesc(shifts, shifts, r);
hold on;
plot(ix, iy, 'kx');

mask_up = createElementMask(Nx_up, Ny_up, dx_up, 2.1e-3, 11.1e-3, ix0_up, iy0_up);
mask_opt = circshift(mask_up, [iy, ix]);

% downsample best mask
mask_corr  = circshift(mask_opt, [floor(upsample_factor/2), floor(upsample_factor/2)]);
mask_down = imresize(mask_corr, sz, 'bilinear', 'Antialiasing', true);
figure;
subplot(2, 2, 1);
hold on;
plot(mask_down(iy0,:), 'k');
plot(ssp_cut(iy0,:)/max(ssp_cut(iy0,:)), 'b');

subplot(2, 2, 2);
hold on;
plot(mask_opt(iy0_up,:) * max(ssp_up( ( (iy0 - 1) * upsample_factor) + 1,: )), 'k');
plot(ssp_up( ( (iy0 - 1) * upsample_factor) + 1,: ), 'b')

subplot(2, 2, 3);
imagesc((mask_down+0.25) .* ssp_cut);
axis image
colormap(getBatlow);

subplot(2, 2, 4);
imagesc((mask_opt+0.25) .* ssp_up);
axis image
colormap(getBatlow);

% Re-assemble the ssp with the border
mask_p_cut = mask_down .* ssp_cut;
mask_2d = zeros(size(ssp));
mask_2d( y_border+1:end-y_border, x_border+1:end-x_border ) = mask_down;


%%
% 
% % sum(ssp, 'all')
% ssp = fracCircShift(ssp, [0, 0.1]);
% % figure;
% % hold on;
% % plot(ssp(30,:));
% % plot(ssp2(30,:));
% % imagesc(ssp2 - ssp);
% % axis image
% 
% % Find optimal mask position
% pt_div      = 10;            % number of divisions per grid point
% delta       = M.dx / pt_div; % the physical length of each shift division
% delta_i     = 1 / pt_div;    % the grid point length of each shift division
% max_shift   = 5*M.dx;          % physical length of largest shift
% Npos        = round(2 * max_shift / delta); % number of shifts
% indexes     = (0:Npos-1) - round((Npos - 1) / 2);
% shifts      = delta_i * indexes;   % the shifts [pts]
% shifts_cart = indexes * delta;     % physical lenght of each shift [m]
% 
% r = zeros(1, Npos);
% % figure;
% % subplot(2, 1, 1);
% % imagesc(ssp);
% % axis image
% % colorbar
% for xdx = 1:Npos
%     for ydx = 1%1:Npos
%         x_shift = shifts(xdx);
%         y_shift = shifts(ydx);
%         mask_shift = fracCircShift(mask, [0, x_shift]);
% 
%         mask_p = mask_shift .* ssp;
%         
% %         subplot(2, 1, 2);
% %         cla(gca);
% %         imagesc(mask_p, [min(ssp(:)), max(ssp(:))]);
% %         axis image
% %         colorbar
% %         drawnow
% %         pause(0.1);
%         r(xdx) = sum(mask_p, 'all');% / sum(mask_shift, 'all');
%     end
% end
% figure;
% plot(shifts, r);
% 
% %%
% 
% 
% [~,i_max] = max(r, [], 'all');
% [y_max, x_max] = ind2sub(size(r), i_max);
% y_shift = shifts(y_max);
% x_shift = shifts(x_max);
% mask = fracCircShift(mask, [y_shift, x_shift]);
% 
% figure;
% imagesc(shifts, shifts, r);
% xlabel('x shift [pts]');
% ylabel('y-shift [pts]');
% hold on
% plot(x_shift, y_shift, 'kx');
% 
% %%
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Preview
% figure;
% subplot(3, 1, 1);
% imagesc(mask);
% axis image
% set(gca, 'YDir', 'normal')
% 
% subplot(3, 1, 2);
% plot(mask(iy0,:))
% title('x cross section')
% 
% subplot(3, 1, 3);
% plot(mask(:,ix0))
% title('y cross section')
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ssp slice before and after masking
% 
% figure;
% plot(ssp(round(M.Ny/2),:));
% 
% % Preview the mask effect on SSP
% mask_p = mask .* ssp;
% 
% % Preview
% figure;
% subplot(3, 1, 1);
% imagesc(mask_p);
% axis image
% set(gca, 'YDir', 'normal')
% title('Unfiltered mask');
% 
% subplot(3, 1, 2);
% plot(mask_p(iy0,:))
% hold on
% plot(ssp(iy0,:))
% % xline(xl, 'k--');
% % xline(xr, 'k--');
% legend({'Masked Pressure', 'Unmasked Pressure'});
% 
% subplot(3, 1, 3);
% plot(mask_p(:,ix0))
% hold on
% plot(ssp(:,ix0))
% % xline(yb, 'k--');
% % xline(yt, 'k--');
% legend({'Masked Pressure', 'Unmasked Pressure'});

%%
% Replicate mask for all time steps
mask_3d = repmat(mask_2d, [1, 1, M.Nt]);

% Apply mask to source plane
mask_source = mask_3d .* M.source_p;
%%
%%
% Project forward to validation plane

% load validation scan data
filename       = 'field_scan_probe_A_channel_8_validation_postprocessed';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
V = load(input_filename, 'pressure', 'z_pos', 'time_axis', 'Nt', 'Nx', 'Ny', 'dt', 'c_water');

% NOTE: this does not need to be rounded to the nearest dx
z_proj = M.source_z + abs(V.z_pos) - abs(M.z_pos);

% % ------------------------------------------------------------------------
% % Show that the effect of grid expansion is only 0.25% max error reduction,
% % at the expense of large memory. GE won't be used.
% as_options = {'GridExpansion', 0, 'Plot', 0, 'Reverse', 0, 'DataCast', 'single', 'DataRecast', 1};
% [~, val_p] = angularSpectrum(mask_source, M.dx, M.dt, z_proj, M.c_water, as_options{:});
% 
% as_options = {'GridExpansion', 200, 'Plot', 0, 'Reverse', 0, 'DataCast', 'single', 'DataRecast', 1};
% [~, val_p_ge] = angularSpectrum(mask_source, M.dx, M.dt, z_proj, M.c_water, as_options{:});
% 
% diff = (val_p - val_p_ge);
% norm_max_err = max(abs(diff), [], 3) / max(val_p(:));
% figure;
% subplot(3, 1, 1);
% imagesc(sum(val_p.^2, 3));
% axis image
% colorbar
% colormap(getBatlow);
% title('Pressure squared integral without grid expansion')
% 
% subplot(3, 1, 2);
% imagesc(sum(val_p_ge.^2, 3));
% axis image
% colorbar
% colormap(getBatlow);
% title('Pressure squared integral with grid expansion')
% 
% subplot(3, 1, 3);
% imagesc(1e2*norm_max_err);
% axis image
% colorbar
% colormap(getBatlow);
% title('Difference [%]')
% % ------------------------------------------------------------------------

%%

as_options = {'GridExpansion', 0, 'Plot', 0, 'Reverse', 0, 'DataCast', 'single', 'DataRecast', 1};
[~, val_p] = angularSpectrum(mask_source, M.dx, M.dt, z_proj, M.c_water, as_options{:});

%%

% Compare measured and validation at the centre frequency 1.2MHz
Fs = 1 / M.dt;
% [a_meas, ~, f_meas] = extractAmpPhase(val_p, Fs, 1.2e6, 'Dim', 3);
[f_meas, as_meas] = spect(val_p, 1/M.dt, 'FFTLength', 4000, 'Dim', 3);
[f_val, as_val] = spect(V.pressure, 1/V.dt, 'FFTLength', 4000, 'Dim', 3);
% [a_val, ~, f_val]   = extractAmpPhase(V.pressure, Fs, 1.2e6, 'Dim', 3);
[f_c, i_c] = findClosest(f_meas, 1.2e6);

%%
ac_meas = as_meas(:,:,i_c);
ac_val = as_val(:,:,i_c);

ac_meas = ac_meas / max(ac_meas(:));
ac_val = ac_val / max(ac_val(:));

figure;
subplot(3, 2, 1);
imagesc(ac_meas);
axis image
colorbar
colormap(getBatlow);
title('Measured, reprojected');

subplot(3, 2, 3);
imagesc(ac_val);
axis image
colorbar
colormap(getBatlow);
title('Validation');

diff = ac_val - ac_meas;
max_norm_err = 1e2* max((diff), [], 3) / max(ac_val(:));

subplot(3, 2, 5);
imagesc(max_norm_err);
axis image
colorbar
colormap(getBatlow);
title('Difference [%]');

subplot(3, 2, 2);
hold on;
plot(ac_val(round(M.Ny/2), :));
plot(ac_meas(round(M.Ny/2), :));
legend({'Validation', 'Measured'});

subplot(3, 2, 4);
hold on;
plot(ac_val(:, round(M.Nx/2)));
plot(ac_meas(:, round(M.Nx/2)));
legend({'Validation', 'Measured'});

% %%
% % Compare measured and validation using sum of squared pressure
% % (normalised)
% ssp_meas = sum(val_p.^2, 3);
% ssp_val  = sum(V.pressure.^2, 3);
% 
% ssp_meas = ssp_meas / max(ssp_meas(:));
% ssp_val = ssp_val / max(ssp_val(:));
% 
% figure;
% subplot(3, 2, 1);
% imagesc(ssp_meas);
% axis image
% colorbar
% colormap(getBatlow);
% title('Measured, reprojected');
% 
% subplot(3, 2, 3);
% imagesc(ssp_val);
% axis image
% colorbar
% colormap(getBatlow);
% title('Validation');
% 
% diff = ssp_val - ssp_meas;
% max_norm_err = 1e2* max((diff), [], 3) / max(ssp_val(:));
% 
% subplot(3, 2, 5);
% imagesc(max_norm_err);
% axis image
% colorbar
% colormap(getBatlow);
% title('Difference [%]');
% 
% subplot(3, 2, 2);
% hold on;
% plot(ssp_val(round(M.Ny/2), :));
% plot(ssp_meas(round(M.Ny/2), :));
% legend({'Validation', 'Measured'});
% 
% subplot(3, 2, 4);
% hold on;
% plot(ssp_val(:, round(M.Nx/2)));
% plot(ssp_meas(:, round(M.Nx/2)));
% legend({'Validation', 'Measured'});

