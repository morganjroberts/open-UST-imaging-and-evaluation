% DESCRIPTION:
%     Backprojecting measurement plane back to the source plane
%
% INPUT DATA FILENAMES:
%     <data-dir>\__________.mat      
%
% EXPERIMENTAL PARAMETERS:
%     <repo-dir>\_____.md
%
% ABOUT:
%     author:      - Morgan Roberts
%     last update: - 24/11/22
% 
% close all
clearvars

[~, data_dir]    = getRepoDataPath();
scan_data_folder = [data_dir, filesep, 'field_scans'];
filename         = 'field_scan_probe_A_all_channels_postprocessed';
input_filename   = [scan_data_folder, filesep, filename, '.mat'];
load(input_filename, 'pressure', 'time_axis', 'Ny', 'Nx', 'c_water', 'x_pos', 'y_pos', 'z_pos', 'dx', 'dt');

plotMeasurementPlane(pressure, x_pos, y_pos)

% Create vector of target plane positions
Nsource_planes = 5;
offset         = -4; % offset (+ve is in the 'reverse' direction, same as scan tank)
z0             = abs(round(z_pos / dx) * dx);
source_planes  = offset + (1:Nsource_planes) - ceil(Nsource_planes/2);
z_vec          = z0 + (source_planes * dx);

% Visually verify all three datasets with quick backpropagation
options      = {'GridExpansion', 0, 'Plot', 0, 'Reverse', 1};
pressure_max = angularSpectrum(pressure, dx, dt, z_vec, c_water, options{:});
pressure_max = flip(pressure_max, 3);
sliceViewer(pressure_max, [dx, dx, dx]);

% -------------------------------------------------------------------------
% Manual selection of best source plane based on image sharpness
source_plane = 3;
source_z     = z_vec(source_plane);

% Backpropagation without grid expansion
options            = {'GridExpansion', 0, 'Plot', 0, 'Reverse', 1};
[~, pressure_time] = angularSpectrum(pressure, dx, dt, source_z, c_water, options{:});

% Backpropagation with grid expansion
options               = {'GridExpansion', Nx, 'Plot', 0, 'Reverse', 1};
[~, pressure_time_ge] = angularSpectrum(pressure, dx, dt, source_z, c_water, options{:});

% Demonstrate that grid expansion is not required
diff_time    = abs(pressure_time_ge - pressure_time);
norm_max_err = 1e2 * max(diff_time, [], 3) / max(pressure_time_ge(:));

figure;
imagesc(norm_max_err);
c = colorbar;
axis image
colormap(getBatlow);
ylabel(c, 'Max difference [%]');
drawnow

% -------------------------------------------------------------------------
% Backpropagate all three datasets in time to the source plane and save
options = {'GridExpansion', 0, 'Plot', 0, 'Reverse', 1};

filename         = 'field_scan_probe_A_all_channels_postprocessed';
input_filename   = [scan_data_folder, filesep, filename, '.mat'];
backpropagate(input_filename, source_z, options);

filename         = 'field_scan_probe_E_all_channels_postprocessed';
input_filename   = [scan_data_folder, filesep, filename, '.mat'];
backpropagate(input_filename, source_z, options);

filename         = 'field_scan_probe_F_all_channels_postprocessed';
input_filename   = [scan_data_folder, filesep, filename, '.mat'];
backpropagate(input_filename, source_z, options);

function backpropagate(input_filename, source_z, options)

load(input_filename, 'pressure', 'c_water', 'dx', 'dt');
[~, source_p]    = angularSpectrum(pressure, dx, dt, source_z, c_water, options{:});
save(input_filename, 'source_p', 'source_z', '-append');

end




