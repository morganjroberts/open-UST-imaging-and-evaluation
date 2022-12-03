close all
clearvars

[~, data_dir]    = getRepoDataPath();
scan_data_folder = [data_dir, filesep, 'field_scans'];

% -------------------------------------------------------------------------
% Load the datasets and backpropagate to a series of planes
filename       = 'field_scan_probe_A_all_channels_postprocessed';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
load(input_filename, 'source_p', 'pressure', 'time_axis', 'Ny', 'Nx', 'Nt', 'source_z', 'c_water', ...
                                    'x_pos', 'y_pos', 'z_pos', 'dx', 'dt');

plotMeasurementPlane(source_p, x_pos, y_pos)
%
[~,Ilin] = max(abs(source_p), [], 'all');
[iy, ix, it] = ind2sub(size(source_p), Ilin);

ix = 98;
truncatePressureInTime(source_p, time_axis, time_axis(end)*1e6, x=ix, y=iy, Clip=0.1);
ylim(gca, [-600,600]);
