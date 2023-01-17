function stats = computeStats(freqs, as, dbThresh)
%COMPUTESTATS calculate mean/std centre-freq, bandwidth, and amplitude.
%
% DESCRIPTION:
%      computeStats calculates the distribution of centre_frequency,
%      fractional bandwidth, and amplitude-at-centre-frequency, for a
%      dataset of amplitude spectra corresponding to Ntdx transmitters.
%      The mean and standard deviation for these quantities are stored. The
%      bandwidth cutoff frequencies are defined as the first time the
%      amplitude spectrum falls below the defined threshold, when moving
%      away from the centre-frequency peak in both directions. These
%      threshold-cutting positions are stored for the mean amplitude
%      spectrum for plotting later.
%
% INPUTS:
%     freqs          - [numeric] frequency vector with Nf frequencies [Hz]
%     as             - [numeric] array with size (Nf, Ntdx) containing the
%                      amplitude spectrum for each transmitter [Pa]
%     dbThresh       - [numeric] threshold in dB for the fractional
%                      bandwidth definition (should be negative)
%
% OUTPUTS:
%     stats          - structure containing the mean and standard deviation
%                      statistics
%
% ABOUT:
%     author         - Morgan Roberts
%     date           - 16th November 2022


% extract centre freq (with maximum amplitude) for each transmitter ------
[~,I]   = max(as, [], 1);
fc      = freqs(I);
mean_fc = mean(fc);
std_fc  = std(fc);

% extract amplitude distribution at mean centre freq----------------------
[fc_cl, I] = findClosest(freqs, mean_fc);
warning(['Actual centre frequency is different from mean_fc by ', num2str(abs(fc_cl-mean_fc)), ' Hz.'])
Afc      = as(I,:);
mean_Afc = mean(Afc);
std_Afc  = std(Afc);

% Compute fractional bandwidth distribution -----------------------------
Ntx      = size(as, 2);
FBW      = zeros(1, Ntx);

for idx = 1:Ntx
    as_spect     = as(:,idx);
    [~,I]        = max(as_spect);
    fc           = freqs(I);
    thresh       = max(as_spect) * 10^(dbThresh/20);

    % Find where spectrum falls below threshold either side of fc (if it rises
    % above threshold again in different frequency range, this is ignored)
    state_change = diff(as_spect > thresh);
    i_low        = find(state_change(1:I), 1, 'last') + 1;
    i_high       = find(state_change(I+1:end), 1, 'first') + I;
    x_low        = freqs(i_low);
    x_high       = freqs(i_high);
    FBW(idx)     = 1e2 * (x_high - x_low) / fc;
end

% repeat for the mean amplitude spectrum to find the points to plot later
as_mean      = mean(as, 2);
[~,I]        = max(as_mean);
thresh       = max(as_mean) * 10^(dbThresh/20);
state_change = diff(as_mean > thresh);
i_low        = find(state_change(1:I), 1, 'last') + 1;
i_high       = find(state_change(I+1:end), 1, 'first') + I;

% Extract coordinates for plotting later
y_low  = as_mean(i_low);
y_high = as_mean(i_high);
mean_y = mean([y_low, y_high]);
y_low = mean_y;
y_high = mean_y;
x_low  = freqs(i_low);
x_high = freqs(i_high);

% compute mean and sttandard deviation fractional bandwith
mean_FBW = mean(FBW);
std_FBW  = std(FBW);

% Store statistics ------------------------------------------------------
stats.dbThresh    = dbThresh;
stats.mean_fc     = mean_fc;
stats.std_fc      = std_fc;
stats.mean_FBW    = mean_FBW;
stats.std_FBW     = std_FBW;
stats.mean_Afc    = mean_Afc;
stats.std_Afc     = std_Afc;
stats.std_Afc_pc  = 1e2 * std_Afc / mean_Afc;
stats.mean.x_low  = x_low;
stats.mean.x_high = x_high;
stats.mean.y_low  = y_low;
stats.mean.y_high = y_high;

end