function [mean_p_trace, mean_p_spect, fig] = plotMeanAndVariation(freqs, as, time_axis, pressure, stats, traceXlim, traceYlim, spectXlim, spectYlim, spectMode)
%PLOTMEANANDVARIATION plots time-pressure and frequency-amplitude
%
% DESCRIPTION:
%     plotMeanAndVariation presents the time-pressure and frequency-amplitude
%     relationship as a subplot. The input dataset is assumed to have
%     information from Ntdx multiple transmitters, so the mean impulse response
%     is shown, in addition to a fill representing the entire range of the
%     measured data.
%
% INPUTS:
%     freqs          - [numeric] frequency vector with Nf frequencies [Hz]
%     as             - [numeric] array with size (Nf, Ntdx) containing the
%                      amplitude spectrum for each transmitter [Pa]
%     time_axis      - [numeric] time vector with Nt samples [s]
%     pressure       - [numeric] pressure array with size (Nt, Ntdx)
%                      containing the hydrophone pressure trace for each
%                      transmitter [Pa]
%     stats          - structure, output from computeStats()
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
%     spectMode      - [char] scale for the amplitude spectrum plots. Set
%                      to 'db' for decibel and 'linear' for linear.
%
% OUTPUTS:
%      mean_p_trace  - [numeric] pressure trace representing the mean
%                      impulse response of all transmitters in the dataset,
%                      has size (1, Nt) [Pa]
%      mean_p_spect  - [numeric] amplitude spectrum representing the mean
%                      amplitude spectrum of all transmitters in the
%                      dataset, has size (1, Nf) [Pa]
%      fig           - figure handle for the plot, to allow saving later
%
% ABOUT:
%      author        - Morgan Roberts
%      date          - 16th November 2022

% Calculate size of vectors
Nt = size(pressure, 1);
Nf = length(freqs);

% storage vectors for 0/0.5/1 quantiles
pressure_quant  = zeros(Nt, 3);
amplitude_quant = zeros(Nf, 3);

mean_p_trace = mean(pressure, 2);
mean_p_spect = mean(as, 2);

% Calculate quantiles (entire range and median)
quant = [0, 0.5, 1];
for fdx = 1:Nf
    amplitude_quant(fdx,:) = quantile(as(fdx,:), quant);
end
for tdx = 1:Nt
    pressure_quant(tdx,:) = quantile(pressure(tdx,:), quant);
end

% calculate the upper/lower bound vectors for the fill plot
f_fill        = [freqs(:)', flip(freqs(:)')];
t_fill        = [time_axis(:)', flip(time_axis(:)')];
pressure_fill = [pressure_quant(:,3)', flip(pressure_quant(:,1)')];
as_fill       = [amplitude_quant(:,3)', flip(amplitude_quant(:,1)')];

fill_col = [175, 238, 238]/255;

% Create figure with specified handle
fig = figure;

% Time domain plot --------------------------------------------------------
subplot(2, 1, 1);
hold on

% plot mean pressure trace
h1 = plot(1e6*time_axis, mean_p_trace*1e-3, 'k-', 'linewidth', 2);
% plot entire range of measured data
h2 = fill(1e6*t_fill, 1e-3*pressure_fill, fill_col);
set(h2, 'EdgeColor','none');
uistack(h2,'bottom');

xlabel('Time [\mus]');
ylabel('Pressure [kPa]');
title('Impulse Response');

set(gca,'Box','on');
set(gca, 'Layer', 'Top')

ylim(traceYlim)
xlim(traceXlim);   

legend([h1, h2], {'Mean', 'Measured Variation'});


% Frequency domain plot ---------------------------------------------------
% if required, change amplitude data to db, instead if linear
if strcmp(spectMode, 'db')
    ref_val           = max(mean_p_spect); % this may need to be an input argument later
    mean_p_spect      = 20 * log10(mean_p_spect / ref_val);
    as_fill           = 20 * log10(as_fill / ref_val);
    stats.mean.y_low  = 20 * log10(stats.mean.y_low / ref_val);
    stats.mean.y_high = 20 * log10(stats.mean.y_high / ref_val);
    as_fill(isinf(as_fill)) = -80; % catch for values of as_fill that are zero, (-Inf in db)
end

subplot(2, 1, 2);
hold on

% plot mean pressure amplitude spectrum
plot(freqs*1e-6, mean_p_spect, 'k-', 'linewidth', 2);

% plot entire range of measured data
h2 = fill(1e-6*f_fill, as_fill, fill_col);
set(h2, 'EdgeColor','none');
uistack(h2,'bottom');

% plot resonance and fractional bandwidth
h5 = xline(stats.mean_fc*1e-6, 'k--');
h6 = plot(1e-6*[stats.mean.x_low, stats.mean.x_high], ...
          [stats.mean.y_low, stats.mean.y_high], ...
          'r-', 'marker', '.', 'markersize', 12);

if strcmp(spectMode, 'db')
    yline(0, 'k--');
end

xlabel('Frequency [MHz]');
ylabel('Pressure [Pa]');
title('Amplitude Spectrum');

set(gca,'Box','on');
set(gca, 'Layer', 'Top');

xlim(spectXlim);
ylim(spectYlim);

legend([h5, h6], {['Mean transducer resonance', newline, ...
    'frequency f_0 = ', num2str(stats.mean_fc*1e-6, 3), ' MHz'], ...
    ['Mean ', num2str(stats.dbThresh), ' dB Fractional Bandwidth: ', num2str(stats.mean_FBW, 2), ' %']}, ...
    'location', 'southeast');

end