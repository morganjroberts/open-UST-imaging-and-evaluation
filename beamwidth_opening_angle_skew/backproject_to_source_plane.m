% FILENAME:
%     backproject_to_source_plane.m
%
% DESCRIPTION:
%     Backprojecting the measurement plane back to the source plane using
%     the angularSpectrum method, for three datasets from probes UST-A,
%     UST_E, and UST_F. First, the exact source plane z-position is
%     estimated by backprojecting to a series of planes surrounding the
%     expected z-position, and visually choosing the plane where the MIP
%     (in time) is sharpest. 
%     Next, it is shown that it is not necessary to
%     use grid expansion for the back-propagation since the errors are
%     ~0.01% relative to the maximum pressure.
%     Finally, the entire pressure time series are backpropagated to the
%     source plane in one step for each dataset, and this 3D volume is
%     appended to the input data file.
%
% INPUT DATA FILENAMES:
%     <data-dir>\field_scans\_______.mat      
%
% EXPERIMENTAL PARAMETERS:
%     <repo-dir>\field_scans\_______.md
%
% ABOUT:
%     author:      - Morgan Roberts
%     last update: - 24/11/22
% 

close all
clearvars

[~, data_dir]    = getRepoDataPath();
scan_data_folder = [data_dir, filesep, 'field_scans'];

% -------------------------------------------------------------------------
% Load the datasets and backpropagate to a series of planes
filename       = 'field_scan_probe_A_all_channels_postprocessed';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
load(input_filename, 'pressure', 'time_axis', 'Ny', 'Nx', 'Nt', 'c_water', ...
                                    'x_pos', 'y_pos', 'z_pos', 'dx', 'dt');
plotMeasurementPlane(pressure, Xpos=x_pos, Ypos=y_pos)

% Create vector of plane positions surrounding the expected z-position
Nsource_planes = 5;
offset         = -4; % offset (+ve is in the 'reverse' direction, same as scan tank)
z0             = abs(round(z_pos / dx) * dx);
source_planes  = offset + (1:Nsource_planes) - ceil(Nsource_planes/2);
z_vec          = z0 + (source_planes * dx);

% Backpropagation to series of source planes using angular spectrum method
options      = {'Plot', 0, 'Reverse', 1};
pressure_max = angularSpectrum(pressure, dx, dt, z_vec, c_water, options{:});
pressure_max = flip(pressure_max, 3); % Required since reverse mode is used
sliceViewer(pressure_max, [dx, dx, dx]);

% Manual selection of best source plane based on MIP image sharpness
source_plane = 3;
source_z     = z_vec(source_plane);

% -------------------------------------------------------------------------
% % Assessment of the effect of grid expansion on this pressure data
% 
% % Backpropagation without grid expansion
% options            = {'Plot', 0, 'Reverse', 1};
% [~, pressure_time] = angularSpectrum(pressure, dx, dt, source_z, c_water, options{:});
% 
% % Backpropagation with grid expansion
% options               = {'GridExpansion', round(Nx/2), 'Plot', 0, 'Reverse', 1};
% [~, pressure_time_ge] = angularSpectrum(pressure, dx, dt, source_z, c_water, options{:});
% 
% % Demonstrate that grid expansion is not required due to very low error
% diff_time    = abs(pressure_time_ge - pressure_time);
% norm_max_err = 1e2 * max(diff_time, [], 3) / max(pressure_time_ge(:));
% 
% figure;
% imagesc(norm_max_err);
% c = colorbar;
% axis image
% colormap(getBatlow);
% ylabel(c, 'Max difference [%]');
% drawnow

% -------------------------------------------------------------------------
% Backpropagate all datasets to the source plane and append to data file
options = {'GridExpansion', 0, 'Plot', 0, 'Reverse', 1};

filename = 'field_scan_probe_A_all_channels_postprocessed';
backpropagate(scan_data_folder, filename, source_z, options);

filename = 'field_scan_probe_E_all_channels_postprocessed';
backpropagate(scan_data_folder, filename, source_z, options);

filename = 'field_scan_probe_F_all_channels_postprocessed';
backpropagate(scan_data_folder, filename, source_z, options);

function backpropagate(data_dir, filename, source_z, as_options)

% Load the pressure measurement plane
input_filename = [data_dir, filesep, filename, '.mat'];
load(input_filename, 'pressure', 'c_water', 'time_axis', 'dx', 'dt', 'Nt');

% Backpropagate to source_z
[~, source_p] = angularSpectrum(pressure, dx, dt, source_z, c_water, as_options{:});

% Centre the source waveform in time and adjust the source time axis
i_shift  = Nt / 2;
source_p = circshift(source_p, [0, 0, i_shift]);
source_t = time_axis - (source_z / c_water) - (i_shift * dt);

% Save the source plane
save(input_filename, 'source_p', 'source_t', 'source_z', '-append');

end




