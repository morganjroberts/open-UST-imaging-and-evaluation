function [time_axis, mean_p_trace, freqs, mean_p_spect] = processImpulseResponse(input_dir, input_filename, cutoff_f, options)
%PROCESSIMPULSERESPONSE Manages the entire transmit-impulse-response processing workflow
%
% DESCRIPTION:
%     processImpulseResponse accepts a single pressure dataset.
%     First the data is filtered to remove the high frequency and DC content.
%     Next, the x-position corresponding to the beam axis is located for
%     everyelement. All other non-beam-axis traces are then discarded.
%     Next, the waveforms are windowed, pre/post padded, and aligned in
%     time. The amplitude spectra are then calculated using an FFT.
%     Finally, everything is plotted, and relevant quantities are
%     calculated.
%
% INPUTS:
%     input_dir      - absolute path to the input data directory
%     input_filename - file containing calibrated pressure traces
%                      (without the .mat suffix)
%     cutoff_f       - [numeric] lowpass cutoff frequency [Hz]
%
% OPTIONAL INPUTS:
%     ExtraPlot      - [boolean] whether to produce extra plots showing
%                      each of the processing stages
%     spectMode      - [char] scale for the amplitude spectrum plots. Set
%                      to 'db' for decibel and 'linear' for linear.
%     dbThresh       - [numeric] threshold in dB for the fractional
%                      bandwidth definition (should be negative)
%     traceXlim      - [char or numeric] set as 'time_axis' to set the
%                      XLim to the start/end of the time axis.
%                      Otherwise, use a 2-element numeric vector [us]
%     traceYlim      - [numeric] 2-element vector setting the pressure
%                      YLim of the trace plot [kPa]
%     spectXlim      - [numeric] 2-element vector setting the frequency
%                      XLim of the spectrum plot. Defaults to the cutoff
%                      frequency [MHz]
%     spectYlim      - [numeric] 2-element vector setting the amplitude
%                      YLim of the spectrum plot [Pa or dB]
%
% OUTPUTS:
%     time_axis      - [numeric] calculated time axis vector corresponding
%                      to the pressure data in the input argument pressure
%                      [s]
%     mean_p_trace   - [numeric] mean pressure signal [Pa]
%     freqs          - [numeric] calculated frequency vector [Hz]
%     mean_p_spect   - [numeric] calculated mean amplitude spectrum [Pa]
%
%
% ABOUT:
%     author         - Morgan Roberts
%     date           - 16th November 2022

arguments
    input_dir
    input_filename
    cutoff_f
    options.ExtraPlot = 0;
    options.spectMode = 'linear';
    options.dbThresh  = -6;
    options.traceXlim = 'time_axis';
    options.traceYlim = [-Inf, Inf];
    options.spectXlim = [0, cutoff_f * 1e-6];
    options.spectYlim = [-Inf, Inf];
end

tic;
disp('Loading pressure data ...');
load([input_dir, input_filename, '.mat'], 'pressure', 'dt', 't0');
disp(['Completed in ', num2str(toc), ' s']);

% create time axis
Nt        = size(pressure, 2);
time_axis = t0 + (0:dt:(Nt - 1) * dt);

% default behaviour for XLim is to cut exactly to the time axis extents
if strcmp(options.traceXlim, 'time_axis')
    options.traceXlim = time_axis([1, end]) * 1e6;
end

% apply pzt polarity correction depending on filename
switch input_filename
    case 'transmit_impulse_response_probe_q'
        pzt_polarity = logical([0,0,1,0,0,0,1,1,0,1,1,1,1,1,1,0]);
    case 'transmit_impulse_response_probes_AEFG'
        load('..\pzt_polarity.mat', 'pzt_polarity');
        pzt_polarity = logical(pzt_polarity([1:16, 65:112]));
    case 'transmit_impulse_response_probe_A22uH'
        load('..\pzt_polarity.mat', 'pzt_polarity');
        pzt_polarity = logical(pzt_polarity(1:16));
end

% correct for polarity differences
pressure(:,:,pzt_polarity) = -pressure(:,:,pzt_polarity);

% filter the input
pressure_filt = filterTraces(pressure, dt, cutoff_f, Plot=options.ExtraPlot);

% locate the x-position that is on the beam axis
pressure_on_axis = locateBeamAxis(pressure_filt, Plot=options.ExtraPlot);

% align the pressure in time
[pressure_shift_pad, pressure_shift] = alignInTime(pressure_on_axis, Plot=options.ExtraPlot);

% compute amplitude spectrum
[freqs, as] = computeSpect(pressure_shift_pad, dt);

% Display stats of centre freq, fractional bandwidth, amplitude at fc
stats = computeStats(freqs, as, options.dbThresh);
disp(input_filename);
disp(stats);

% plot pressure-time and amplitude-frequency with mean and variation
[mean_p_trace, mean_p_spect, fig] = ... 
    plotMeanAndVariation(freqs, as, time_axis, pressure_shift, stats, ...
    options.traceXlim, options.traceYlim, options.spectXlim, ...
    options.spectYlim, options.spectMode);

% Title and save figure
title(gca, input_filename)
set(fig,'renderer','Painters');
figure_filename = [pwd, '\figures\', input_filename];
print(fig, figure_filename, '-depsc2');
print(fig, figure_filename, '-dsvg');
end