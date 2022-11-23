% FILENAME:
%      transmit_receive_directional_response_analysis.m
%
% DESCRIPTION:
%     Code for analysing the watershot receive data from the open-UST
%     transducer ring array. Assessing the effect of emission and incidence
%     angle on the mean and standard deviation of the received amplitude
%     (at centre freq). Tx/Rx rays are placed into bins based on the angle
%     of emission/incidence relative to the element normals. Within each
%     bin, signals are extracted and amplitude spectra are calculated.
%     Statistics are then computed within each bin and plotted.
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

% Sampling frequency for FFT below
Fs = settings.acqSamplingFrequency;

% Load predicted element positions
load([repo_dir, filesep, 'calibration', filesep, 'ideal_element_positions.mat'], 'angle_TxRx_to_TxNorm', 'angle_RxTx_to_RxNorm');

% Set angular bin width, discretise the emission/incidence angle matrices
step       = 5;
thetaR_bin = round(angle_RxTx_to_RxNorm / step) * step;
thetaT_bin = round(angle_TxRx_to_TxNorm / step) * step;

% Use mask for restricting to the angles that we have receive data for
thetaR_bin(~mask) = NaN;
thetaT_bin(~mask) = NaN;

% Extract the emission/incidence angles corresponding to bin centres.
bins = unique(thetaR_bin(~isnan(thetaR_bin)));
Nbin = length(bins);

% ------------------------------------------------------------------------
% Plot of the emission/incidence angles of the ray relative to the
% transmitting/receiving element normal
figure;
subplot(1, 2, 1);
h = imagesc(thetaT_bin);
axis image;
c = colorbar;
colormap(getBatlow());
title('Emission');
xlabel('Receivers');
ylabel('Transmitters');
ylabel(c, 'Angle [deg]');
set(h, 'AlphaData', ~isnan(thetaT_bin));

subplot(1, 2, 2);
h = imagesc(thetaR_bin);
axis image;
c = colorbar;
colormap(getBatlow());
title('Incidence');
xlabel('Receivers');
ylabel('Transmitters');
ylabel(c, 'Angle [deg]');
set(h, 'AlphaData', ~isnan(thetaR_bin));
% ------------------------------------------------------------------------

% Loop over each unique emission/incidence angular bin combination
% - extract the grouped-and-aligned receive data for that bin
% - calculate the amplitude spectra
% - compute statistics

% Storage vectors for statistics
Nel      = NaN * ones(Nbin);
mean_Afc = NaN * ones(Nbin);
std_Afc  = NaN * ones(Nbin);
mean_fc  = NaN * ones(Nbin);
std_fc   = NaN * ones(Nbin);

for btx = 1:Nbin
    for brx = 1:Nbin
        % Create emission/incidence angle mask
        t_mask = abs(thetaT_bin - bins(btx)) < 1e-6;
        r_mask = abs(thetaR_bin - bins(brx)) < 1e-6;
        tr_mask = and(t_mask, r_mask);

        % Check if there is data avaiable for this angular bin combination
        Nel(btx, brx) = sum(tr_mask(:));
        if Nel(btx, brx) == 0
            continue
        end

        % Extract the receive data traces in this angular bin combination
        g_data = groupSimilarSignals(fa_data, tr_mask, Align=true, RestrictXcorr=6, ExtraPlot=false);

        % Calculate the amplitude spectra for each receive data trace
        [freqs, as] = spect(g_data, Fs, 'Dim', 2, 'FFTLength', size(g_data, 2) * 6);
        
        % Compute and store statistics for the current angular bin
        dbThresh           = -10;
        stats              = computeStats(freqs, as', dbThresh);
        mean_Afc(btx, brx) = stats.mean_Afc;
        std_Afc(btx, brx)  = stats.std_Afc;
        mean_fc(btx, brx)  = stats.mean_fc;
        std_fc(btx, brx)   = stats.std_fc;        
    end
end

figure; 
subplot(1, 3, 1);
h = imagesc(Nel);
set(h, 'AlphaData', ~isnan(std_Afc));
colorbar
title('N');
subplot(1, 3, 2);
h = imagesc(mean_fc);
set(h, 'AlphaData', ~isnan(std_Afc));
colorbar
title('Mean f_c');
subplot(1, 3, 3);
h = imagesc(std_fc);
set(h, 'AlphaData', ~isnan(std_Afc));
colorbar
title('Std f_c');

% ------------------------------------------------------------------------
% Plot the effect of emission and incidence angle on the mean and standard
% deviation of the recieve amplitude at the centre frequency.
fig1 = figure;
subplot(1, 2, 1);
h = imagesc(bins, bins, 20*log10(mean_Afc/max(mean_Afc(:))));
axis image;
c = colorbar;
colormap(getBatlow());
ylabel('Emission [deg]');
xlabel('Incidence [deg]');
set(h, 'AlphaData', ~isnan(mean_Afc));
ylabel(c, 'Mean Amplitude [dB]')
xlim(bins([1, end]));
ylim(bins([1, end]));

subplot(1, 2, 2);
h = imagesc(bins, bins, 1e2*std_Afc/max(mean_Afc(:)));
axis image;
c = colorbar;
colormap(getBatlow());
ylabel('Emission [deg]');
xlabel('Incidence [deg]');
set(h, 'AlphaData', ~isnan(std_Afc));
ylabel(c, 'Standard Deviation [%]')
xlim(bins([1, end]));
ylim(bins([1, end]));

set(fig1, 'Position', [680 817 817 281]);

set(fig1,'renderer','Painters');
figure_filename = [repo_dir, filesep, 'transmit_receive_response', filesep, 'figures', filesep, 'tx_rx_directional_response_mean_std'];
print(fig1, figure_filename, '-depsc2');
print(fig1, figure_filename, '-dsvg');
