% FILENAME:
%     beam_width_analysis.m
%
% DESCRIPTION:
%     Script to compute the elevation beam width of N=48 UST elements.
%
% INPUT DATA FILENAMES:
%     <data-dir>\field_scans\_______.mat      
%
% ABOUT:
%     author:      - Morgan Roberts
%     last update: - 6/12/22

close all
clearvars

[~, data_dir]    = getRepoDataPath();
scan_data_folder = [data_dir, filesep, 'field_scans'];

tic;
disp('Loading data ...');
filename       = 'isolated_reprojected_beam_data_probes_AEF';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
load(input_filename, 'beam_data', 'freqs', 'c_water', 'dx', 'Ny', 'Nx', 'Fs', 'dt', 'z_proj', 'z_targ', 'Nfreqs');
disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');
%%
% Calculate x and y position vectors
x_pos = 0:dx:(Nx-1) * dx;
x_pos = x_pos - x_pos( round(Nx / 2) );
y_pos = 0:dx:(Ny-1) * dx;
y_pos = y_pos - y_pos( round(Ny / 2) );

% Setup loop parameters
Nel = size(beam_data, 5);
zdx = 5;

% Setup storage
ele_resp    = zeros(Ny, Nfreqs, Nel);
ele_resp_fc = zeros(Ny, Nel);
ycs         = zeros(1, Nel);
beamwidth   = zeros(1, Nel);

% Threshold for defining the edge of the beam
thresh = -6;

% Loop over each element
for edx = 1:Nel
    % Extract the data plane at z=110mm for the current element
    % volume_data has dims (Ny, Nx, Nfreqs)
    volume_data = squeeze( beam_data(:, :, :, zdx, edx) );

    % Find the centroid and edges of the -6dB bounding box in 3D
    [X, ~, F] = locateThresholdRegion(volume_data, thresh);

    % Extract and store the elevational response spectrum at the centroid
    % This has dims (Ny, Nf)
    ele_resp(:, :, edx) = squeeze( volume_data(:, X.ic, :) );

    % Extract and store the elevational cross section at f_c
    % This has dims (1, Ny)
    ele_xsec = volume_data(:, X.ic, F.ic);
    ele_resp_fc(:, edx) = ele_xsec;

    % Compute the beamwidth
    beamwidth(edx) = fwhm(ele_xsec, dx, 0);
  
end

% % Correct for the beam skew by shifting 
% figure;
% plot(y_pos, squeeze(ele_lines(100,:,:)) )
% 
% Npad = 15;
% pad = NaN * ones(Nfreqs, Npad, Nel);
% ele_lines_pad = cat(2, pad, ele_lines, pad);
% ele_lines_shift = NaN * ones(size(ele_lines_pad));
% i_centre = round(Ny / 2);
% for edx = 1:Nel
%     i_shift = i_centre - ycs(edx);
%     ele_lines_shift(:,:,edx) = circshift(ele_lines_pad(:,:,edx), [0, i_shift]);
% end
% 
% % Remove the NaN padding
% ele_lines_shift = ele_lines_shift(:, Npad+1:end-Npad, :);
% 
% figure;
% plot( squeeze(ele_lines_shift(100,:,:)) )
% 
% ele_lines = ele_lines_shift;

% Calculate the mean elevational response and the standard deviation
mean_bw = mean(beamwidth);
std_bw  = std(beamwidth);
disp(['Elevation Beamwidth mean = ', num2str(1e3*mean_bw, 3),  ' mm, std = ', num2str(1e3*std_bw, 3), ' mm']);

% Calculate the mean elevational response, convert to dB and compute
% standard deviation
mean_resp = mean(ele_resp, 3);
max_val   = max(mean_resp(:));
mean_dB   = 20 * log10(mean_resp / max_val);
std_resp  = std(ele_resp, 0, 3);
std_pc    = 1e2 * std_resp / max_val;

% ------------------------------------------------------------------------
% Plot the mean elevational response and the standard deviation

% Find the -6dB edges of the mean elevational response for plotting
[~, Y, ~] = locateThresholdRegion( permute( mean_resp, [1, 2, 3] ), thresh );
plot_thresh = -30;

fig = figure;
subplot(1, 2, 1);
hold on;
imagesc(y_pos*1e3, freqs*1e-6, mean_dB', [plot_thresh, 0]);
c = colorbar;
set(gca, 'YDir', 'normal');
colormap(getBatlow);
axis square
ylabel('Frequency [MHz]');
xlabel('y-position [mm]');
ylabel(c, 'Pressure [dB]');
% xline(y_pos(Y.imin)*1e3, 'k--', 'linewidth', 1.5);
% xline(y_pos(Y.imax)*1e3, 'k--', 'linewidth', 1.5);
xlim( 1e3 * y_pos([1, end]) );
ylim( 1e-6 * freqs([1, end]) );

subplot(1, 2, 2);
imagesc(y_pos*1e3, freqs*1e-6, std_pc');
c = colorbar;
set(gca, 'YDir', 'normal');
colormap(getBatlow);
axis square
ylabel('Frequency [MHz]');
xlabel('y-position [mm]');
ylabel(c, 'Standard Deviation [%]');
xlim( 1e3 * y_pos([1, end]) );
ylim( 1e-6 * freqs([1, end]) );

set(fig, 'Position', [404 278 932 602]);
set(fig,'renderer','Painters');
filename = [pwd, '\figures\elevational_response'];
print(fig, filename, '-dsvg');
print(fig, filename, '-depsc2');

% ------------------------------------------------------------------------
% Plot the beamwidth histogram

face_colour = [175, 238, 238]/255;

figure;
bin_edges = 15:0.25:17.5;
h2 = histogram(beamwidth*1e3, bin_edges);
set(h2,'facecolor',face_colour);
xlabel('Beamwidth [mm]');
ylabel('Counts');
axis square
xlim(bin_edges([2, end]))

% set(fig,'renderer','Painters');
% filename = [pwd, '\figures\beamwidth_histogram'];
% print(fig, filename, '-dsvg');
% print(fig, filename, '-depsc2');

