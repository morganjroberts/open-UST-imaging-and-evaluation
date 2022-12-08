% FILENAME:
%     skew_analysis.m
%
% DESCRIPTION:
%     Script to compute the elevation and lateral beam axis skews of N=48
%     elements. First, the bulk module skew is estimated from the field
%     scan data. Next, 5 planes are loaded for each of the elements,
%     containing pressure amplitude spectra over an XY grid of points.
%     Within each plane the centroid is located. A line is fit in 3D to the
%     5 centroids, defining the beam axis. The angle of this line is
%     calculated in the elevation and lateral planes and stored. The bulk
%     lateral skew of the modules (due to scan plane misalignment) is then
%     removed.
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

% Estimate the bulk skew from each module and allocate to elements
filename      = 'field_scan_probe_A_all_channels_postprocessed';
thetaA        = calculateBulkLateralSkew(scan_data_folder, filename, ExtraPlot=true);
filename      = 'field_scan_probe_E_all_channels_postprocessed';
thetaE        = calculateBulkLateralSkew(scan_data_folder, filename, ExtraPlot=true);
filename      = 'field_scan_probe_F_all_channels_postprocessed';
thetaF        = calculateBulkLateralSkew(scan_data_folder, filename, ExtraPlot=true);
bulk_lat_skew = [thetaA*ones(1, 16), thetaE*ones(1, 16), thetaF*ones(1, 16)];

tic;
disp('Loading data ...');
filename       = 'isolated_reprojected_beam_data_probes_AEF';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
load(input_filename, 'beam_data', 'freqs', 'c_water', 'dx', 'Ny', 'Nx', 'Fs', 'dt', 'z_proj', 'z_targ', 'Nfreqs');
disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');

% Calculate x and y position vectors
x_pos = 0:dx:(Nx-1) * dx;
x_pos = x_pos - x_pos( round(Nx / 2) );
y_pos = 0:dx:(Ny-1) * dx;
y_pos = y_pos - y_pos( round(Ny / 2) );

% Setup loop parameters
Nel = size(beam_data, 5);
Nz  = size(beam_data, 4);

% Setup storage
beam_x   = zeros(Nz, Nel);
beam_y   = zeros(Nz, Nel);
skew_ele = zeros(1, Nel);
skew_lat = zeros(1, Nel);

% Threshold for defining the edge of the beam
thresh = -6;

% Replicate z-coordinates across elements
beam_z = repmat( z_targ(:), [1, Nel] );

% Loop over each element
for edx = 1:Nel
    for zdx = 1:Nz
        % Extract the data plane
        % volume_data has dims (Ny, Nx, Nfreqs)
        volume_data = squeeze( beam_data(:, :, :, zdx, edx) );
    
        % Find the centroid of the -6dB bounding box in 3D
        [X, Y, ~] = locateThresholdRegion(volume_data, thresh, Round=false);

        % Store the centroid coordinates
        beam_x(zdx, edx) = x_pos( floor(X.ic) ) + (dx * abs( X.ic - floor(X.ic) ));
        beam_y(zdx, edx) = y_pos( floor(Y.ic) ) + (dx * abs( Y.ic - floor(Y.ic) ));
    end

    % Fit 3D line https://www.codefull.net/2015/06/3d-line-fitting/
    points     = [squeeze( beam_x(:,edx) ), ...
                  squeeze( beam_y(:,edx) ), ..., 
                  squeeze( beam_z(:,edx) )];
    p0         = mean(points, 1);                % Find mean of points
    subtracted = bsxfun(@minus, points, p0);     % Subtract mean
    [~, ~, V]  = svd(subtracted);                % Perform SVD
    dir        = V(:, 1);                        % Find direction vector

    % Extract parametric coefficients
    a  = dir(1);
    b  = dir(2);
    c  = dir(3);

    % Use 3D parametric fitting coefficients to define beam axis angles   
    skew_lat(edx) = atand(a / c); 
    skew_ele(edx) = atand(b / c);  

end

% Calculate absolute skews, with bulk lateral skew corrected
abs_skew_lat = abs(skew_lat - bulk_lat_skew);
abs_skew_ele = abs(skew_ele);

% Calculate means and standard deviations
skew_ele_mean = mean(abs_skew_ele);
skew_ele_std  = std(abs_skew_ele);
skew_lat_mean = mean(abs_skew_lat);
skew_lat_std  = std(abs_skew_lat);

disp(['Elevation skew mean = ', num2str(skew_ele_mean, 3),  ' deg, std = ', num2str(skew_ele_std, 3), ' deg']);
disp(['Lateral skew mean = ', num2str(skew_lat_mean, 4),  ' deg, std = ', num2str(skew_lat_std, 3), ' deg']);

% -----------------------------------------------------------------------
% Plot the skews

face_colour = [175, 238, 238]/255;

figure;
subplot(1, 2, 1);
h = histogram(abs_skew_ele, 10);
set(h,'facecolor',face_colour);
xlabel('Skew [deg]');
ylabel('Counts');
xlim([0, Inf]);
title('Elevation Skew');
axis square

subplot(1, 2, 2);
h = histogram(abs_skew_lat, 10);
set(h,'facecolor',face_colour);
xlabel('Skew [deg]');
ylabel('Counts');
xlim([0, Inf]);
title('Lateral Skew');
axis square
