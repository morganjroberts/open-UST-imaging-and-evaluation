% FILENAME:
%     post_process_validation_scans.m
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

filename         = 'line_scanX_probe_A_channel_8_validation';
postprocess(scan_data_folder, filename)

filename         = 'line_scanY_probe_A_channel_8_validation';
postprocess(scan_data_folder, filename)

filename         = 'field_scan_probe_A_channel_8_validation';
postprocess(scan_data_folder, filename)

% ------------------------------------------------------------------------
function postprocess(scan_data_folder, filename)
    input_filename   = [scan_data_folder, filesep, filename, '.mat'];
    load(input_filename, 'pressure', 'time_axis', 'Ny', 'Nx', 'dt', 'Nt');

    % Make linescan into a fieldscan
    switch filename
        case {'line_scanX_probe_A_channel_8_validation'}
            pressure = permute(pressure, [3, 1, 2]); 
            [pressure, time_axis, Nt] = truncatePressureInTime(pressure, time_axis, 97.5, Clip=0.02, TukeyParam=0.05, Plot=true);
        case  'line_scanY_probe_A_channel_8_validation'
            pressure = permute(pressure, [1, 3, 2]); 
            [pressure, time_axis, Nt] = truncatePressureInTime(pressure, time_axis, 97.5, Clip=0.02, TukeyParam=0.05, Plot=true);
           
    end
    disp(size(pressure));
    cutoff_f  = 2e6;   % set to 2 MHz since the spatial sampling step doesn't support higher than ~2.1MHz

    % Low pass filter the pressure
    pressure = applyFilterVolume(pressure, dt, cutoff_f, 3, RemovePad=true, PadLength=2, ExtraPlot=true);

    % Post-pad the signal and adjust the time axis to make the same size as the
    % other datasets
    pad_factor = 1;
    pad        = zeros(Ny, Nx, Nt*pad_factor);
    pressure   = cat(3, pressure, pad);
    time_axis  = [time_axis, (time_axis(end) + dt + (0:dt:(Nt - 1) * dt))]; % only valid for postpad = 1!
    Nt         = length(time_axis);
    
    % Save the post-processed data to the same directory with a new filename
    output_filename = [scan_data_folder, filesep, filename, '_postprocessed.mat'];
    if ~exist(output_filename, 'file')
        copyfile(input_filename, output_filename)
    end
    save(output_filename, 'pressure', 'time_axis', 'Nt', '-append');

end