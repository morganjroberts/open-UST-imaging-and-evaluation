% FILENAME:
%      transmit_receive_on_axis_response_analysis.m
%
% DESCRIPTION:
%     Code for analysing the watershot receive data from the open-UST
%     transducer ring array. Traces from the Tx/Rxs that are exactly
%     opposite one another are grouped, and their amplitude spectra are
%     calculated. The centre frequency, bandwidth, and amplitude-at-centre
%     freq are then computed. The results are plotted to help assess the
%     effect of interelement variation.
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

% Calculate time and frequency parameters
Nt = size(fa_data, 3);
Fs = settings.acqSamplingFrequency;
dt = 1/Fs;
time_axis = 0:dt:((Nt - 1) * dt);

% Load mask showing which elements are opposite one another
load([repo_dir, filesep, 'calibration', filesep, 'ideal_element_positions.mat'], 'is_opposite');

% Plot the mask
figure;
imagesc(is_opposite);
axis image
title('Mask for opposite elements');

% Extract the aligned data for all rays where tx/rx are exactly opposite
g_data = groupSimilarSignals(fa_data, is_opposite, Align=true, RestrictXcorr=6, ExtraPlot=false);

% Calculate amplitude spectra for these traces
[freqs, as] = spect(g_data, Fs, 'Dim', 2, 'FFTLength', Nt * 6);

% Compute statistics on the amplitude spectra
dbThresh = -12;
stats = computeStats(freqs, as', dbThresh);
disp(stats);
disp(stats.mean);

% Plot the traces and spectra with mean and variation shown
i_start = find(sum(g_data, 1), 1, 'first') - 10;
i_end   = find(sum(g_data, 1), 1, 'last') + 10;
[mean_p_trace, mean_p_spect, fig] = plotMeanAndVariation(freqs, as', time_axis, g_data', stats, 1e6*time_axis([i_start, i_end]), [-Inf,Inf], [0.4, 2.5], [-40, 2], 'db');

% Save the figure
set(fig,'renderer','Painters');
figure_filename = [repo_dir, filesep, 'transmit_receive_response', filesep, 'figures', filesep, 'tx_rx_on_axis_response'];
print(fig, figure_filename, '-depsc2');
print(fig, figure_filename, '-dsvg');
