% FILENAME:
%      transmit_receive_response_analysis.m
%
% DESCRIPTION:
%     previous version that uses different angular binning approach
%
% NOTE:
%     This script doesn't use the raw receive data - this should have
%     previously been modified by extract_first_arrival_wave_packet.m
%
% INPUT DATA FILENAMES:
%     <data-dir>\UST_dataset1\First_arrival_capture_UST_watershot_2-11-22 16-18-28.mat Modified UST watershot dataset
%
% FIGURE OUTPUT DIRECTORY:
%     <repo-dir>\transmit_receive_response\figures
%
% ABOUT:
%     author:      - Morgan Roberts
%     last update: - 15/11/22

close all
clearvars

% Load the watershot data
[repo_dir, data_dir] = getRepoDataPath;
dataset  = [data_dir, filesep, 'UST_dataset1'];
filename = 'First_arrival_capture_UST_watershot_2-11-22 16-18-28.mat';
tic;
disp('Loading watershot data ...');
load([dataset, filesep, filename], 'fa_data', 'mask', 'settings');
disp(['Completed in ', num2str(toc),  's']);
fprintf('\n');

Fs = settings.acqSamplingFrequency;

% Load predicted element positions
load([repo_dir, filesep, 'calibration', filesep, 'ideal_element_positions.mat'], 'element_positions', 'distance_matrix', 'angle_TxRx_to_TxNorm', 'angle_TxNorm_to_RxNorm', 'is_opposite');
%%
load('lateral_directional_response.mat')
figure;
plot(theta, response);

angle_RxTx_to_RxNorm = angle_TxRx_to_TxNorm';

dir_resp = NaN * ones(size(distance_matrix));
Ntdx = size(fa_data, 1);
Nrdx = size(fa_data, 2);
Nt = size(fa_data, 3);
for tdx = 1:Ntdx
    for rdx = 1:Nrdx
        if mask(tdx, rdx)
            tx_ang = angle_TxRx_to_TxNorm(tdx, rdx);
            rx_ang = angle_RxTx_to_RxNorm(tdx, rdx);
            tx_dir_resp = interp1(theta, response, tx_ang);
            rx_dir_resp = interp1(theta, response, rx_ang);
            dir_resp(tdx, rdx) = tx_dir_resp * rx_dir_resp;
        end
    end
end
figure;
subplot(1, 2, 1); imagesc(dir_resp); axis image; colorbar
subplot(1, 2, 2); plot(squeeze(dir_resp(1,:)));

%%
step = 0.05;
dir_binned = round(dir_resp / step) * step;
% fig1 = figure;
imagesc(dir_binned); axis image; colorbar; drawnow
bins = flip(0:step:1);
Nbin = length(bins);

binned_traces  = NaN * ones(Nbin, 3000);
binned_spectra = NaN * ones(Nbin, 8000);
A_at_fc = NaN * ones(Nbin, 1);
std_Afc = NaN * ones(Nbin, 1);
mean_Afc = NaN * ones(Nbin, 1);

for bdx = 1:Nbin
    disp(['dir_val ', num2str(bins(bdx))]);
    dir_mask = abs(dir_binned - bins(bdx)) < 1e-6;
    g_data = groupSimilarSignals(fa_data, dir_mask, Normalise=false, RestrictXcorr=6, ExtraPlot=false);
    if length(g_data) == 1 && g_data == 0
        continue
    end
    [freqs, as] = spect(g_data, Fs, 'Dim', 2, 'FFTLength', Nt * 6);
    
    i_start = find(sum(g_data, 1), 1, 'first') - 10;
    i_end   = find(sum(g_data, 1), 1, 'last') + 10;

    Nf = length(freqs);

    dbThresh = -10;
    stats = computeStats(freqs, as', dbThresh);
    time_axis = 1e-6*(1:Nt);
    close gcf
    [mean_p_trace, mean_p_spect, fig] = plotMeanAndVariation(freqs, as', time_axis, g_data', stats, [i_start, i_end], [-20,20], [0, 2.5], [0, 200], 'linear');
    binned_traces(bdx, 1: (i_end-i_start+1) ) = mean_p_trace(i_start:i_end);
    binned_spectra(bdx, 1:Nf) = mean_p_spect;
    std_Afc(bdx) = stats.std_Afc;
    mean_Afc(bdx) = stats.mean_Afc;
%     drawnow
%     a = input('enter');
end

% Cut arrays to size (excess space preallocated earlier)
i_stop = find(~sum(~isnan(binned_traces), 1), 1) - 1;
binned_traces = binned_traces(:,1:i_stop);
i_stop = find(~sum(~isnan(binned_spectra), 1), 1) - 1;
binned_spectra = binned_spectra(:,1:i_stop);

figure;
hold on;
valid_bins = find(~isnan(binned_traces(:,1)));
B = length(valid_bins);
for bdx = 1:B
    subplot(1, 4, 2);
    hold on;
    rx_data = binned_traces(valid_bins(bdx),:);
    rx_data = 0.5* rx_data / max(abs(binned_traces(:)));
    plot(-rx_data + bdx, 'k');
    xlim([1, 160]);
    xlabel('Time [samples]');
    set(gca, 'YDir', 'reverse')
    ylim([0, B+0.5]);
    ylabel('Bin Index');

    subplot(1, 4, 3);
    hold on;
    rx_data = binned_spectra(valid_bins(bdx),:);
    rx_data = rx_data / max(abs(binned_spectra(:)));
    plot(freqs*1e-6, -rx_data + bdx, 'k');    
    xlim([0, 3])
    ylim([0, B+0.5]);
    xlabel('Frequency [MHz]');
    set(gca, 'YDir', 'reverse')

end
% subplot(1, 4, 3)
% imagesc(freqs*1e-6,1:B,binned_spectra(valid_bins,:));
% xlim([0, 3])
% xlabel('Frequency [MHz]');

subplot(1, 4, 1)
h=plot((bins(valid_bins)), 1:B, 'k');
h.Marker = '.';
ylim([0, B+0.5]);
set(gca, 'YDir', 'reverse')
xlabel('Directional Response [au]');

subplot(1, 4, 4)
hold on
plot(mean_Afc(valid_bins), 1:B, 'k');
plot(mean_Afc(valid_bins) - std_Afc(valid_bins), 1:B, 'k--');
plot(mean_Afc(valid_bins) + std_Afc(valid_bins), 1:B, 'k--');
ylim([0, B+0.5]);
xlim([0, Inf])
set(gca, 'YDir', 'reverse')
xlabel('Amplitude at fc');
