% DESCRIPTION:
%     Post-processing field scans for 
%
% INPUT DATA FILENAMES:
%     <data-dir>\__________.mat        module Q (without acoustic matching layers)
%
% EXPERIMENTAL PARAMETERS:
%     <repo-dir>\_____.md
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
% Preview effect of post-processing parameters

t_cut       = 58;
border      = 4;
tukey_param = 0.05;
cutoff_f    = 2e6; % set to 2 MHz since the spatial sampling step doesn't support higher than ~2.1MHz

plotMeasurementPlane(pressure, x_pos, y_pos)

% Inspect traces and cut hydrophone reflection
ix = round(Nx/2);
iy = round(Ny/2);
[pressure, time_axis, Nt] = removeTipReflection(pressure, time_axis, t_cut, Clip=0.02, TukeyParam=tukey_param, x=ix, y=iy, Plot=true);

% Window the measurement plane in space
pressure = spatialWindowMeasurementPlane(pressure, border);

plotMeasurementPlane(pressure, x_pos, y_pos)

% ------------------------------------------------------------------------
% Perform post processing for all files
filename = 'field_scan_probe_A_all_channels';
postprocessFieldScan(scan_data_folder, filename, t_cut, border, cutoff_f, TukeyParam=tukey_param)

filename = 'field_scan_probe_E_all_channels';
postprocessFieldScan(scan_data_folder, filename, t_cut, border, cutoff_f, TukeyParam=tukey_param)

filename = 'field_scan_probe_F_all_channels';
postprocessFieldScan(scan_data_folder, filename, t_cut, border, cutoff_f, TukeyParam=tukey_param)


% DO THIS FOR ALL THREE AT ONCE? plane must be in same place
% stage 1 backproject using angular spectrum using careful settings - max p
% only - select the best source plane
% stage 2, backproject entire time series to just the source plane
% save the pressure volumes
% stage 3 locate and isolate (window) the sources. Save
% stage 4 re-project to 5(?) planes and save
% stage 5 compute beamwidth, opening angle, skew, save
% stage 6 visualise stats

function postprocessFieldScan(input_dir, input_filename, t_cut, border, cutoff_f, options)

arguments
    input_dir
    input_filename
    t_cut
    border
    cutoff_f
    options.TukeyParam
end

% load data
tic;
disp('Loading pressure data ...');
file_path = [input_dir, filesep, input_filename, '.mat'];
load(file_path, 'pressure', 'time_axis', 'dt');
disp(['Completed in ', num2str(toc), ' s']);

% Remove hydrophone tip reflection
[pressure, time_axis, Nt] = removeTipReflection(pressure, time_axis, t_cut, Clip=0.02, TukeyParam=options.TukeyParam, Plot=true);

% Apply a spatial window to the measurement plane
pressure = spatialWindowMeasurementPlane(pressure, border);

% Low pass filter the measurement data
pressure = applyFilterVolume(pressure, dt, cutoff_f, 3, RemovePad=true, PadLength=2, ExtraPlot=true);

% Save the data to the same directory with a new filename
output_filename = [input_dir, filesep, input_filename, '_postprocessed.mat'];
if ~exist(output_filename, 'file')
    copyfile(file_path, output_filename)
end
save(output_filename, 'pressure', 'time_axis', 'Nt', '-append');

end

function [pressure_win, time_axis_cut, Nt] = removeTipReflection(pressure, time_axis, t_end, options)

arguments
    pressure
    time_axis
    t_end
    options.Plot       = true;
    options.x          = round(size(pressure, 2)/2);
    options.y          = round(size(pressure, 1)/2);
    options.Clip       = 0.02;
    options.TukeyParam = 0.05;
end

tic;
disp('Removing Hydrophone Tip Reflection...');

ix = options.x;
iy = options.y;

Nx = size(pressure, 2);
Ny = size(pressure, 1);

slice_y = squeeze(pressure(iy, :, :));
trace   = squeeze(pressure(iy, ix, :));

lw = 2;

% Find exact discretised cutting position
[t_cut, Nt] = findClosest(time_axis, t_end*1e-6);

% Trim the pressure data and update other parameters
pressure_cut  = pressure(:,:,1:Nt);
time_axis_cut = time_axis(1:Nt);

% Create a right-sided sided window to smooth the cut, replicate in 3D
win                            = getWin(Nt, 'Tukey', 'Param', options.TukeyParam);
win( 1:find(win == 1, 1) - 1 ) = 1;
win                            = permute(win, [2, 3, 1]);
win                            = repmat(win, [Ny, Nx, 1]);

% Apply the window
pressure_win = pressure_cut .* win;

if options.Plot
    figure;
    subplot(3, 1, 1);
    imagesc(time_axis*1e6, 1:Nx, slice_y, [-options.Clip, options.Clip] * max(slice_y(:)));
    xlim(time_axis([1, end]) * 1e6);
    xlabel('Time [us]');
    ylabel('x-position [samples]');
    hold on;
    xline(t_cut*1e6, 'b', 'linewidth', lw);
    yline(ix, 'k', 'linewidth', lw);
    colormap(getBatlow);
    
    slice_y2 = squeeze(pressure_win(iy, :, :));
    trace2   = squeeze(pressure_win(iy, ix, :));
    
    subplot(3, 1, 2);
    imagesc(time_axis_cut*1e6, 1:Nx, slice_y2, [-options.Clip, options.Clip] * max(slice_y2(:)));
    xlim(time_axis([1, end]) * 1e6);
    xlabel('Time [us]');
    ylabel('x-position [samples]');
    hold on;
    xline(t_cut*1e6, 'b', 'linewidth', lw);
    yline(ix, 'k', 'linewidth', lw);
    colormap(getBatlow);
    
    subplot(3, 1, 3);
    hold on
    plot(time_axis*1e6, trace*1e-3, 'k');
    plot(time_axis_cut*1e6, trace2*1e-3, 'r');
    xlim(time_axis([1, end]) * 1e6);
    xlabel('Time [us]');
    ylabel('Pressure [kPa]');
    hold on;
    xline(t_cut*1e6, 'b--', 'linewidth', lw);
    legend({'Original Trace', 'Trimmed Trace'})

    drawnow
end

disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');

end

function pressure_win = spatialWindowMeasurementPlane(pressure, border)

tic;
disp('Apply spatial window to measurement plane ...');

Ny = size(pressure, 1);
Nx = size(pressure, 2);
Nt = size(pressure, 3);

% Create window and extract the ring up and ring down regions
win = getWin(border*2+1, 'Tukey', 'Param', 1);
ru  = win(1:border);
rd  = win( (border + 2):end );

% Create new windows, controlling the width of the middle region
win_row = [ru; ones(Ny - (2 * border), 1); rd];
win_col = [ru; ones(Nx - (2 * border), 1); rd];

% Use the outer product to create the 2D window
win_2d = win_row * win_col';

% Replicate the window for all time steps
win_3d = repmat(win_2d, [1, 1, Nt]);

% Apply the window
pressure_win = pressure .* win_3d;

disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');

end

