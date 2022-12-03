% FILENAME:
%     post_process_field_scans.m
%
% DESCRIPTION:
%     Post-processing field scans to prepare them for back-propagation.
%     First, the data is truncated to remove the reflection between the
%     transducer face and hydrophone tip. Windowing is used to smooth this
%     transition. Windowing is also applied in both spatial dimensions at
%     every time step to prepare for the FFTs during backpropagation.
%     Next, the time series are filtered to remove high frequency
%     content not supported by the spatial grid step size. Finally, the
%     time series are post-padded. To save the data, the input data file is
%     copied, and the pressure, time_axis, and Nt variables are overwritten
%     with their post-processed versions. All other variables are kept the
%     same.
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

close all
clearvars

[~, data_dir]    = getRepoDataPath();
scan_data_folder = [data_dir, filesep, 'field_scans'];
filename         = 'field_scan_probe_A_all_channels';
input_filename   = [scan_data_folder, filesep, filename, '.mat'];
load(input_filename, 'pressure', 'time_axis', 'Ny', 'Nx', 'x_pos', 'y_pos');

% ------------------------------------------------------------------------
% First preview the effect of the post-processing parameters

t_cut       = 58;
border      = 4;
tukey_param = 0.05;
cutoff_f    = 2e6; % set to 2 MHz since the spatial sampling step doesn't support higher than ~2.1MHz

% Visualise measurement plane before processing
plotMeasurementPlane(pressure, x_pos, y_pos)

% Inspect traces, truncate pressure array to remove hydrophone tip reflection
ix = round(Nx/2);
iy = round(Ny/2);
[pressure, time_axis, Nt] = truncatePressureInTime(pressure, time_axis, t_cut, Clip=0.02, TukeyParam=tukey_param, x=ix, y=iy, Plot=true);

% Window the measurement plane in space
pressure = spatialWindowMeasurementPlane(pressure, border);

% Visualise measurement plane after processing
plotMeasurementPlane(pressure, x_pos, y_pos)

% ------------------------------------------------------------------------
% Now, perform the post processing for all files

filename = 'field_scan_probe_A_all_channels';
postprocessFieldScan(scan_data_folder, filename, t_cut, border, cutoff_f, TukeyParam=tukey_param, PostPadFactor=1)

filename = 'field_scan_probe_E_all_channels';
postprocessFieldScan(scan_data_folder, filename, t_cut, border, cutoff_f, TukeyParam=tukey_param, PostPadFactor=1)

filename = 'field_scan_probe_F_all_channels';
postprocessFieldScan(scan_data_folder, filename, t_cut, border, cutoff_f, TukeyParam=tukey_param, PostPadFactor=1)
