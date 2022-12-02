function pressure_filt = applyFilterVolume(pressure, dt, filter_param, dim, options)
%APPLYFILTERVOLUME Lowpass-filters the pressure signals with zero-phase
%
% DESCRIPTION:
%      applyFilterVolume takes a 3D set of pressure traces from a hydrophone
%      measurement and applies a zero-phase low pass filter with -80 dB
%      stop band attenuation and a transition width of 0.001 based on the
%      temporal sampling frequency. The DC component is also removed.
%      Padding and windowing is used to prevent the zero-phase filter from
%      producing discontinuities. If required, the unfiltered/filtered
%      signals can be plotted. Nt: number of time samples.
%
% INPUTS:
%     pressure      - [numeric] input pressure array, must have dims: 3 [Pa]
%     dt            - [numeric] time step [s]
%     filter_param  - [numeric] scalar value for the low pass cutoff
%                     frequency [Hz]
%     dim           - [numeric] integer index of the dimension
%                     corresponding to time in the pressure array
%
% OPTIONAL INPUTS:
%     ExtraPlot     - [boolean] Whether to plot an example of a filtered
%                     signal and a comparison of its unfiltered/filtered
%                     amplitude spectrum
%     RemovePad     - [boolean] Symmetrical pre and post padding is added
%                     before filtering. Set this to true to
%                     remove the padding after filtering, before returning
%                     the pressure array
%     PadLength     - [numeric] length of the pre and post padding regions
%                     as an integer multiple of Nt. The padded signal
%                     length is (2 * PadLength + 1) * Nt
%
% OUTPUTS:
%     pressure_filt - [numeric] filtered pressure array [Pa]
%
% ABOUT:
%     author        - Morgan Roberts
%     date          - 16th November 2022 

arguments
    pressure
    dt
    filter_param
    dim
    options.ExtraPlot      = 0;
    options.RemovePad = 1;
    options.PadLength = 2;
end

tic;
disp('Filtering pressure data...');

if length(size(pressure)) ~= 3
    error('Pressure input must have dims: 3');
end

% Work out which array dimensions are not in the time direction
dims          = 1:3;
non_time_dims = dims(dims ~= dim);
dimA          = non_time_dims(1);
dimB          = non_time_dims(2);

% Compute sizes of arrray
Nt = size(pressure, dim);
Na = size(pressure, dimA);
Nb = size(pressure, dimB);

% Permute array to have size (Nt, Na, Nb)
dim_order = [dim, dimA, dimB];
pressure  = permute(pressure, dim_order);

% Initalise storage array with matching dimension order
if options.RemovePad
    output_time_length = Nt;
else
    output_time_length = Nt * ((2 * options.PadLength) + 1);
end
pressure_filt = zeros(output_time_length, Na, Nb);

% Remove DC component of pressure
pressure = removeDCOffset(pressure, 1);

% Set the filter options
Fs        = 1 / dt; 
filt_opts = {'LowPass', 'StopBandAtten', 80, 'TransitionWidth', 0.001, 'ZeroPhase', 1};

% Create a window
win = getWin(Nt, 'Tukey', 'Param', 0.1);

% Create a paddng vector
pad = zeros(Nt * options.PadLength, 1);

% Loop through and Lowpass filter each trace
for adx = 1:Na
    for bdx = 1:Nb
        % Extract trace
        trace = squeeze(pressure(:, adx, bdx));
    
        % Apply window
        trace = trace(:) .* win(:);
    
        % Pre/post pad to prevent zero-phase filter creating discontinuities
        trace_pad  = [pad; trace; pad];
        filt_trace = applyFilter(trace_pad, Fs, filter_param, filt_opts{:});

        % remove the pre/post padding, and window again
        if options.RemovePad
            filt_trace = filt_trace( (Nt * options.PadLength) + 1 : end - (Nt * options.PadLength) );
            filt_trace = filt_trace' .* win;
        end

        % Save the filtered trace
        pressure_filt(:, adx, bdx) = filt_trace;

        % Plot an example
        if options.ExtraPlot && adx == round(Na/2) && bdx == round(Nb/2)
            applyFilter(trace_pad, Fs, filter_param, filt_opts{:}, 'Plot', 1);
            xlim([0, 10]);
            drawnow

            figure;
            hold on;
            if options.RemovePad
                plot(trace, 'k')
            else
                plot(trace_pad, 'k');
            end
            plot(filt_trace, 'r')
            xlabel('Time [samples]');
            ylabel('Pressure [Pa]');
            drawnow
        
        end
    end
end

% Permute back to original dimension order
[~, unsort_dims] = sort(dim_order);
pressure_filt    = permute(pressure_filt, unsort_dims);

disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');

end