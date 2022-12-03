function postprocessFieldScan(input_dir, input_filename, t_cut, border, cutoff_f, options)
%POSTPROCESSFIELDSCAN trunctaes, windows in 3D and filters a pressure measurement plane
%
% DESCRIPTION:
%     postprocessFieldScan loads a pressure measurement plane dataset and
%     post processes it before to prepare it for back-propagation. First,
%     the data is truncated to remove the reflection between the transducer
%     face and hydrophone tip. Windowing is used to smooth this transition.
%     Windowing is also applied in both spatial dimensions at every time
%     step to prepare for the FFTs during backpropagation. Next, the time
%     series are filtered to remove high frequency content not supported by
%     the spatial grid step size. Finally, the time series are post-padded.
%     To save the data, the input data file is copied, and the pressure,
%     time_axis, and Nt variables are overwritten with their post-processed
%     versions. All other variables are kept the same.
%
% USAGE:
%     postprocessFieldScan(input_dir, input_filename, t_cut, border,
%     cutoff_f)
%
% INPUTS:
%     input_dir      - [char] absolute path to the input data directory
%                      (no file separator)
%     input_filename - file containing calibrated pressure traces
%                      (without the .mat suffix)
%     t_cut          - [numeric] the new end time of the pressure data [us]
%     border         - [numeric] integer number of samples defining the
%                      width of the border region
%     cutoff_f       - [numeric] lowpass cutoff frequency [Hz]
%
% OPTIONAL INPUTS:
%     TukeyParam     - [numeric] scalar value between 0 and 1 defining the
%                      taper width of the tukey window used to smooth the
%                      truncation cut.
%     PostPadFactor  - [integer] multiple-of-Nt number of zeros to
%                      post-pad the time traces with after filtering
%
% ABOUT:
%     author         - Morgan Roberts
%     date           - 2/12/22

arguments
    input_dir
    input_filename
    t_cut
    border
    cutoff_f
    options.TukeyParam    = 0.05;
    options.PostPadFactor = 1;
end

% Load data
tic;
disp('Loading pressure data ...');
file_path = [input_dir, filesep, input_filename, '.mat'];
load(file_path, 'pressure', 'time_axis', 'dt', 'Nx', 'Ny');
disp(['Completed in ', num2str(toc), ' s']);

% Remove hydrophone tip reflection
[pressure, time_axis, Nt] = truncatePressureInTime(pressure, time_axis, t_cut, Clip=0.02, TukeyParam=options.TukeyParam, Plot=true);

% Apply a spatial window to the measurement plane
pressure = spatialWindowMeasurementPlane(pressure, border);

% Low pass filter the measurement data
pressure = applyFilterVolume(pressure, dt, cutoff_f, 3, RemovePad=true, PadLength=2, ExtraPlot=true);

% Post-pad the signal and adjust the time axis
pad       = zeros(Ny, Nx, Nt*options.PostPadFactor);
pressure  = cat(3, pressure, pad);
time_axis = [time_axis, (time_axis(end) + dt + (0:dt:(Nt - 1) * dt))];
Nt        = length(time_axis);

% Save the post-processed data to the same directory with a new filename
output_filename = [input_dir, filesep, input_filename, '_postprocessed.mat'];
if ~exist(output_filename, 'file')
    copyfile(file_path, output_filename)
end
save(output_filename, 'pressure', 'time_axis', 'Nt', '-append');

end