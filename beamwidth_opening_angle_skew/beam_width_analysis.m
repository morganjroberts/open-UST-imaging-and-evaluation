close all
clearvars

[~, data_dir]    = getRepoDataPath();
scan_data_folder = [data_dir, filesep, 'field_scans'];


tic;
disp('Loading data ...');
filename       = 'isolated_reprojected_beam_data_ch1';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
load(input_filename, 'beam_data', 'freqs', 'c_water', 'dx', 'Ny', 'Nx', 'Fs', 'dt', 'z_proj', 'z_targ', 'Nfreqs');
disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');

% Extract the plane at z = 110mm, the middle of the ring array
mid_plane = beam_data(:,:,:,end);

% Calculate x and y position vectors
x_pos = 0:dx:(Nx-1) * dx;
x_pos = x_pos - x_pos( round(Nx / 2) );
y_pos = 0:dx:(Ny-1) * dx;
y_pos = y_pos - y_pos( round(Ny / 2) );

%%
% Extract the plane for the largest frequency component
[~,I]        = max(mid_plane, [], 'all');
[~, ~, i_f]  = ind2sub(size(mid_plane), I);
mid_plane_fc = mid_plane(:,:,i_f);

% Threshold the plane [dB rel to max]
mid_plane_dB = 20 * log10(mid_plane_fc / max(mid_plane_fc(:)));
thresh = -6;
mid_plane_thresh = mid_plane_dB;
mid_plane_thresh(mid_plane_thresh < thresh) = NaN;
mask = ~isnan(mid_plane_thresh);

% Compute the location of the weighted centroid of the thresholded plane
props = regionprops(mask, mid_plane_thresh, 'WeightedCentroid');
i_x   = round(props.WeightedCentroid(1));
i_y   = round(props.WeightedCentroid(2));

figure;
h = imagesc(mid_plane_thresh);
set(h, 'AlphaData', mask);
colorbar
colormap(getBatlow);
axis image
xline(i_x, 'k--');
yline(i_y, 'k--');

% ------------------------------------------------------------------------
figure;
subplot(3, 1, 1);
h = imagesc(x_pos*1e3, y_pos*1e3, mid_plane_thresh);
set(h, 'AlphaData', mask);
c = colorbar;
colormap(getBatlow);
axis image
xline(x_pos(i_x)*1e3, 'k--');
yline(y_pos(i_y)*1e3, 'k--');
xlabel('x-position [mm]');
xlabel('y-position [mm]');
ylabel(c, 'Pressure [dB]');
title('z=110mm at f=1.2MHz');

subplot(3, 1, 2);
hold on;
plot(x_pos*1e3, squeeze( mid_plane_dB(i_y,:) ), 'k' );
title('Lateral Plane');
xlabel('x-position [mm]');
ylabel('Pressure [dB]');

subplot(3, 1, 3);
hold on;
plot(y_pos, squeeze( mid_plane_dB(:,i_x) ), 'k' );
title('Elevation Plane');
xlabel('y-position [mm]');
ylabel('Pressure [dB]');

%%
% Use the centroids to extract the elevational and lateral planes

ele_beam = squeeze(mid_plane(:,i_x,:))';
lat_beam = squeeze(mid_plane(i_y,:,:))';

% ele_beam = ele_beam.^2;
% lat_beam = lat_beam.^2;

% Heuristic smoothing for now until mean data is available across all 
k = 7;
ele_beam = movmean(ele_beam, k, 1);
lat_beam = movmean(lat_beam, k, 1);

ele_beam_dB = 20 * log10(ele_beam / max(ele_beam(:)));
lat_beam_dB = 20 * log10(lat_beam / max(lat_beam(:)));

ele_beam_mask = ele_beam_dB > thresh;
lat_beam_mask = lat_beam_dB > thresh;

ele_y_sum = sum(ele_beam_mask, 1);
ele_L     = find(ele_y_sum, 1, 'first');
ele_R     = find(ele_y_sum, 1, 'last');
ele_f_sum = sum(ele_beam_mask, 2);
ele_T     = find(ele_f_sum, 1, 'first');
ele_B     = find(ele_f_sum, 1, 'last');

lat_y_sum = sum(lat_beam_mask, 1);
lat_L     = find(lat_y_sum, 1, 'first');
lat_R     = find(lat_y_sum, 1, 'last');
lat_f_sum = sum(lat_beam_mask, 2);
lat_T     = find(lat_f_sum, 1, 'first');
lat_B     = find(lat_f_sum, 1, 'last');

figure;
subplot(1, 2, 1);
h = imagesc(y_pos*1e3, freqs*1e-6, ele_beam_dB, [-20, 0]);
% set(h, 'AlphaData', ele_beam_mask);
% contourf(y_pos*1e3, freqs*1e-6, ele_beam_dB, -10:1:0);
set(gca, 'YDir', 'normal');
ylabel('Frequency [MHz]');
xlabel('y-position [mm]');
title('Elevation Plane (Squared - TxRx)');
c = colorbar;
colormap(getBatlow);
ylabel(c, 'Pressure [dB]');
xline( y_pos(ele_L)*1e3, 'k--' );
xline( y_pos(ele_R)*1e3, 'k--' );
yline( freqs(ele_T)*1e-6, 'k--' );
yline( freqs(ele_B)*1e-6, 'k--' );
axis square

subplot(1, 2, 2);
h = imagesc(x_pos*1e3, freqs*1e-6, lat_beam_dB, [-20, 0]);
% set(h, 'AlphaData', lat_beam_mask);
% contourf(x_pos*1e3, freqs*1e-6, lat_beam_dB, -10:1:0);
set(gca, 'YDir', 'normal');
ylabel('Frequency [MHz]');
xlabel('x-position [mm]');
title('Lateral Plane');
c = colorbar;
colormap(getBatlow);
ylabel(c, 'Pressure [dB]');
xline( x_pos(lat_L)*1e3, 'k--' );
xline( x_pos(lat_R)*1e3, 'k--' );
yline( freqs(lat_T)*1e-6, 'k--' );
yline( freqs(lat_B)*1e-6, 'k--' );
axis square



