
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
load([scan_data_folder, filesep, filename, '.mat']);

% Low-pass filter and window the measured traces
pressure = applyFilterVolume(pressure, dt, 5e6, 3, RemovePad=true, PadLength=2, ExtraPlot=true);

%%
% data = squeeze( pressure(round(Ny/2), round(Nx/2), :) );
% figure;
% plot(  data);
% spect(data, Fs, 'Plot', 1, 'FFTLength', Nt*10);


% ----------------------------------------------------------------------
% plot slice and trace
% Look at energy coverage of measurement plane
ssp = sum(pressure.^2, 3);
figure;
imagesc(x_pos, y_pos, ssp);
axis image
colormap(getBatlow)
figure;
contourf(ssp/max(ssp(:)), 0:0.1:1);
colorbar;
colormap(getBatlow)
axis image

% zero-phase filter and window traces in time
% window in space
% decide whether cutting is necessary

% DO THIS FOR ALL THREE AT ONCE? plane must be in same place
% stage 1 backproject using angular spectrum using careful settings - max p
% only - select the best source plane

% stage 2, backproject entire time series to just the source plane
% save the pressure volumes

% stage 3 locate and isolate (window) the sources. Save

% stage 4 re-project to 5(?) planes and save

% stage 5 compute beamwidth, opening angle, skew, save

% stage 6 visualise stats


function pressure_filt = filterFieldScanTraces(pressure, dt, filter_param, options)


arguments
    pressure
    dt
    filter_param
    options.ExtraPlot = 0;
end

tic;
disp('Filtering pressure data...');

% Compute sizes of arrray
Ny            = size(pressure, 1);
Nx            = size(pressure, 2);
Nt            = size(pressure, 3);
pressure_filt = zeros(Ny, Nx, Nt);
Fs            = 1 / dt;

% Remove DC component of pressure
pressure = removeDCOffset(pressure, 2);

% Set the filter options
filt_opts = {'LowPass', 'StopBandAtten', 80, 'TransitionWidth', 0.001, 'ZeroPhase', 1};

% Create a window
win = getWin(Nt, 'Tukey', 'Param', 0.03);

% Create a paddng vector
pad = zeros(Nt, 1);

% plot an example if required using builtin plotting from applyFilter
if options.ExtraPlot
    trace = pressure(round(Ny/2),round(Nx/2),:);
    
    % window
    trace = trace(:) .* win(:);

    % pre/post pad to prevent zero-phase filter creating discontinuities
    trace = [pad; trace; pad];

    % apply the filter
    filt_trace = applyFilter(trace, Fs, filter_param, filt_opts{:}, 'Plot', 1);

    % remove the pre/post padding
    filt_trace = filt_trace(Nt+1:end-Nt);
    trace      = trace(Nt+1:end-Nt);

    % window again
    filt_trace = filt_trace .* win';

    xlim([0, 10]);
    figure;
    hold on;
    plot(trace, 'k')
    plot(filt_trace, 'r')
    xlabel('Time [samples]');
    ylabel('Pressure [Pa]');
    drawnow
end

% Loop through and Lowpass filter each trace
for idx = 1:Ny
    for jdx = 1:Nx
        trace = squeeze(pressure(idx,jdx,:));
    
        % window
        trace = trace(:) .* win(:);
    
        % pre/post pad to prevent zero-phase filter creating discontinuities
        trace_pad = [pad; trace; pad];
        filt_trace = applyFilter(trace_pad, Fs, filter_param, filt_opts{:});

        % remove the pre/post padding
        filt_trace = filt_trace(Nt+1:end-Nt);

        % window again
        filt_trace = filt_trace .* win';

        pressure_filt(idx,jdx,:) = filt_trace;
    end
end
disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');


end