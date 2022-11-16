function pressure_filt = filterTraces(pressure, dt, filter_param, options)
%FILTERTRACES Lowpass-filters the pressure signals with zero-phase
%
% DESCRIPTION:
%      filterTraces takes a set of pressure traces from a hydrophone
%      measurement and applies a zero-phase low pass filter with -80 dB
%      stop band attenuation and a transition width of 0.001 based on the
%      temporal sampling frequency. The DC component is also removed.
%      Padding and windowing is used to prevent the zero-phase filter from
%      producing discontinuities. If required, the unfiltered/filtered
%      signals can be plotted. Npos: number of linescan measurement
%      positions. Nt: number of time samples. Ntdx: number of transmitters
%      in the dataset.
%
% INPUTS:
%     pressure      - [numeric] input pressure array with size Npos x Nt x
%                     Ntdx. [Pa]
%     dt            - [numeric] time step [s]
%     filter_param  - [numeric] scalar value for the low pass cutoff
%                     frequency [Hz]
%
% OPTIONAL INPUTS:
%     Plot          - [boolean] Whether to plot an example of a filtered
%                     signal and a comparison of its unfiltered/filtered
%                     amplitude spectrum
%
% OUTPUTS:
%     pressure_filt - [numeric] filtered pressure array with size
%                     (Npos, Nt*3, Ntdx) [Pa]
%
% ABOUT:
%     author        - Morgan Roberts
%     date          - 16th November 2022 

arguments
    pressure
    dt
    filter_param
    options.Plot = 0;
end

tic;
disp('Filtering pressure data...');

% Compute sizes of arrray
Npos          = size(pressure, 1);
Nt            = size(pressure, 2);
Ntdx          = size(pressure, 3);
pressure_filt = zeros(Npos, Nt*3, Ntdx);
Fs            = 1 / dt;

% Remove DC component of pressure
pressure = removeDCOffset(pressure, 2);

% Set the filter options
filt_opts = {'LowPass', 'StopBandAtten', 80, 'TransitionWidth', 0.001, 'ZeroPhase', 1};

% Create a window
win = getWin(Nt, 'Tukey', 'Param', 0.1);

% Create a paddng vector
pad = zeros(Nt, 1);

% plot an example if required using builtin plotting from applyFilter
if options.Plot
    trace = pressure(round(Npos/2),:,1);
    
    % window
    trace = trace(:) .* win(:);

    % pre/post pad to prevent zero-phase filter creating discontinuities
    trace = [pad; trace; pad];

    % apply the filter
    filt_trace = applyFilter(trace, Fs, filter_param, filt_opts{:}, 'Plot', 1);

    % remove the pre/post padding
%     filt_trace = filt_trace(Nt+1:end-Nt);
%     trace      = trace(Nt+1:end-Nt);

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
for idx = 1:Npos
    for jdx = 1:Ntdx
        trace = pressure(idx,:,jdx);
    
        % window
        trace = trace(:) .* win(:);
    
        % pre/post pad to prevent zero-phase filter creating discontinuities
        trace_pad = [pad; trace; pad];
        filt_trace = applyFilter(trace_pad, Fs, filter_param, filt_opts{:});

        % remove the pre/post padding
%         filt_trace = filt_trace(Nt+1:end-Nt);

        pressure_filt(idx,:,jdx) = filt_trace;
    end
end
disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');

end