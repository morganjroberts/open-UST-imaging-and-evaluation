% FILENAME:
%     isolation_reprojection_mask_size_validation.m
%
% DESCRIPTION:
%     It is required to assess the beam shapes of N=48 elements. Instead of
%     perofrming field scans for multiple elements, entire modules were
%     scanned with all elements transmitting, to save time and memory. To
%     estimate the beam shapes of individual elements from a scan with
%     multiple elements being driven, the source pressure distribution of
%     each element must be isolated at the source plane, using a mask.
%     Since the mask design affects the field shape, the mask design must
%     be validated. As a comparison, a field scan was performed where only
%     a single element is transmitting (probe A channel 8). Also in this
%     script, it is shown that grid expansion can be set to zero for the
%     angular spectrum reprojections, to save time and memory, with
%     negligible change in error.
%
% APPROACH:
%     Pre-requisite: Project measurement plane back to source plane
%     (completed in backproject_to_source_plane.m)
%     1. Create an array of different mask sizes. 
%     2. Use the function extractSourceElement.m to extract the source
%     pressure distribution of channel 8 only, for each of the mask sizes.
%     3. For each separate extracted pressure, re-project to the validation
%     plane (z = 80.15mm)
%     4. Calculate the normalised amplitude at 1.2MHz for this re-projected plane
%     5. Also calculate the normalised amplitude at 1.2MHz for the validation scan
%     6. Store the sum of absolute differences between (4) and (5).
%     7. Repeat for all mask sizes. The mask sizes producing the lowest
%     error is optimal, and will be used for the isolation/reprojection
%     pipeline to estimate beamshapes of individual elements.
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

% -------------------------------------------------------------------------
% Load the datasets and backpropagate to a series of planes
filename       = 'field_scan_probe_A_all_channels_postprocessed';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
% M structure is the meausred pressure
M = load(input_filename, 'source_p', 'pressure', 'time_axis', 'Ny', 'Nx', 'Nt', 'source_z', 'c_water', ...
                                    'x_pos', 'y_pos', 'z_pos', 'dx', 'dt');

% load validation scan data
filename       = 'field_scan_probe_A_channel_8_validation_postprocessed';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
V = load(input_filename, 'pressure', 'z_pos', 'time_axis', 'Nt', 'Nx', 'Ny', 'dt', 'c_water');

% For UST probes in UCL scan tank, channel 1 has most +ve x-coordinate
% Visual identification of element centre (ch8) after border cutting
plotMeasurementPlane(M.source_p);
ix0 = 148; 
iy0 = 31; 

% Distance between source plane and validation plane
% NOTE: this does not need to be rounded to the nearest dx
z_proj = M.source_z + abs(V.z_pos) - abs(M.z_pos);

% ------------------------------------------------------------------------
% Show that the effect of grid expansion is only 0.25% max error reduction,
% at the expense of large memory. GE won't be used.
demonstrate_ge = 0;
if demonstrate_ge

    as_options = {'GridExpansion', 0, 'Plot', 0, 'Reverse', 0, 'DataCast', 'single', 'DataRecast', 1};
    [~, val_p] = angularSpectrum(extract_p, M.dx, M.dt, z_proj, M.c_water, as_options{:});
    
    as_options = {'GridExpansion', 200, 'Plot', 0, 'Reverse', 0, 'DataCast', 'single', 'DataRecast', 1};
    [~, val_p_ge] = angularSpectrum(extract_p, M.dx, M.dt, z_proj, M.c_water, as_options{:});
    
    diff = (val_p - val_p_ge);
    norm_max_err = max(abs(diff), [], 3) / max(val_p(:));
    figure;
    subplot(3, 1, 1);
    imagesc(sum(val_p.^2, 3));
    axis image
    colorbar
    colormap(getBatlow);
    title('Pressure squared integral without grid expansion')
    
    subplot(3, 1, 2);
    imagesc(sum(val_p_ge.^2, 3));
    axis image
    colorbar
    colormap(getBatlow);
    title('Pressure squared integral with grid expansion')
    
    subplot(3, 1, 3);
    imagesc(1e2*norm_max_err);
    axis image
    colorbar
    colormap(getBatlow);
    title('Difference [%]')
end
% ------------------------------------------------------------------------

% Setup angular spectrum parameters
as_options   = {'GridExpansion', 0, 'Plot', 0, 'Reverse', 0, 'DataCast', 'single', 'DataRecast', 1};

% Create range of mask sizes, initialise storage
phys_size    = [0.8e-3, 9.8e-3];
offsets      = 0.9e-3:0.05e-3:1.9e-3;
mask_sizes   = [1e-3, 10e-3] + offsets(:);
No           = length(offsets);
sum_norm_err = zeros(1, No);

% Compute amplitude distribution at centre frequency for validation plane
f_c             = 1.2e6;
Fs              = 1 / M.dt;
NFFT            = 4e3;
[f_val, as_val] = spect(V.pressure, 1/V.dt, 'FFTLength', NFFT, 'Dim', 3);
[~, i_c]        = findClosest(f_val, f_c);
ac_val          = as_val(:, :, i_c);
ac_val          = ac_val / max( ac_val(:) );

% Loop over different mask sizes and compute amplitude distribution,
% compare with validation plane, store sum of normalised absolute errors
for odx = 1:No
    mask_size = mask_sizes(odx, :);
    mask_size = [2.35e-3, 11.35e-3];
    
    % Extract the source pressure from a single channel
    extract_p = extractSourceElement(M.source_p, M.dx, [iy0, ix0], ...
        Border=[10,80], UpSample=7, ExtraPlot=false, ...
        PhysicalSize=phys_size, MaskSize=mask_size);
    
    % Project the source to the measurement plane
    [~, val_p] = angularSpectrum(extract_p, M.dx, M.dt, z_proj, M.c_water, as_options{:});
    
    % Compute error between measured and validation at the centre frequency 1.2MHz
    [~, as_meas] = spect(val_p, 1/M.dt, 'FFTLength', NFFT, 'Dim', 3);
    
    % Extract amplitude at desired frequency
    ac_meas = as_meas(:,:,i_c);

    % Normalise to help agreement
    ac_meas = ac_meas / max(ac_meas(:));

    % Compute and store the sum of the normalised errors
    abs_diff          = abs(ac_val - ac_meas);
    sum_norm_err(odx) = sum( abs_diff / max(ac_val(:)), 'all' );

end

fig = figure;
h = plot(offsets*1e3, sum_norm_err, 'k');
h.Marker = '.';
xlabel('Offsets [mm]');
ylabel('Sum of normalised errors');
set(fig, 'renderer', 'Painters');
figure_filename = [pwd, '\figures\', 'isolation_mask_validation_errors'];
savefig(fig, figure_filename);


% Calculate optimal mask design that gives the lowest error
% RESULT: mask_size_opt = [0.00235, 0.01135]
[~,I]         = min(sum_norm_err);
offset_opt    = offsets(I);
mask_size_opt = mask_sizes(I, :);

% Extract the source pressure from a single channel with optimal mask
extract_p = extractSourceElement(M.source_p, M.dx, [iy0, ix0], ...
    Border=[10,80], UpSample=7, ExtraPlot=true, ...
    PhysicalSize=phys_size, MaskSize=mask_size_opt);

% Project the source to the measurement plane
[~, val_p] = angularSpectrum(extract_p, M.dx, M.dt, z_proj, M.c_water, as_options{:});

% Compute error between measured and validation at the centre frequency 1.2MHz
[~, as_meas] = spect(val_p, 1/M.dt, 'FFTLength', NFFT, 'Dim', 3);

% Extract amplitude at desired frequency
ac_meas = as_meas(:,:,i_c);

% Normalise to help agreement
ac_meas = ac_meas / max(ac_meas(:));

% Compute and store the sum of the normalised errors
diff     = ac_val - ac_meas;
norm_err = 1e2 * diff / max(ac_val(:));

%%
% Plot the agreement between measured/reprojected and validation
fig = figure;
subplot(3, 2, 1);
imagesc(ac_meas);
axis image
colorbar
colormap(getBatlow);
title('Measured, reprojected');
xlabel('x-position [samples]');
ylabel('y-position [samples]');

subplot(3, 2, 3);
imagesc(ac_val);
axis image
colorbar
colormap(getBatlow);
title('Validation');
xlabel('x-position [samples]');
ylabel('y-position [samples]');

subplot(3, 2, 5);
imagesc(norm_err);
axis image
c = colorbar;
colormap(getBatlow);
title('Difference [%]');
xlabel('x-position [samples]');
ylabel('y-position [samples]');
ylabel(c, 'Normalised error relative to max value [%]');

mip_val_x  = max(ac_val, [], 1);
mip_meas_x = max(ac_meas, [], 1);
mip_val_y  = max(ac_val, [], 2);
mip_meas_y = max(ac_meas, [], 2);

subplot(3, 2, 2);
hold on;
plot( mip_val_x );
plot( mip_meas_x );
legend({'Validation', 'Measured'});
xlabel('x-position [samples]');
ylabel('Normalised amplitude at f_c');
title('Lateral Plane MIP');

subplot(3, 2, 4);
hold on;
plot( mip_val_y );
plot( mip_meas_y );
legend({'Validation', 'Measured'});
xlabel('y-position [samples]');
ylabel('Normalised amplitude at f_c');
title('Elevation Plane MIP')

set(fig, 'renderer', 'Painters');
figure_filename = [pwd, '\figures\', 'isolation_mask_mips'];
savefig(fig, figure_filename);
