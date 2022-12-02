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
offset         = 2; % offset in the reverse direction
z0             = abs(round(z_pos / dx) * dx);
source_planes  = offset + (0:(Nsource_planes - 1));
z_vec          = z0 + (source_planes * dx);

options = {'GridExpansion', 100, 'Plot', 1, 'Reverse', 1};
pressure_max = angularSpectrum(pressure, dx, dt, z_vec, c_water, options{:});