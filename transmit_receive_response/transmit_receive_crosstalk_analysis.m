% FILENAME:
%     transmit_receive_crosstalk_analysis.m
% 
% DESCRIPTON:

% INPUT DATA FILENAMES:
%     <data-dir>\UST_dataset1\UST_acquisition_UST-watershot1_80V_1avg_475tgc_1cyc_1.4045MHz_02-Nov-2022 16-18-28.mat  UST watershot dataset
%
% OUTPUT DATA FILENAMES:
%     <data-dir>\UST_dataset1\First_arrival_capture_UST_watershot_2-11-22 16-18-28.mat                                Modified UST watershot dataset
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
filename = 'UST_acquisition_UST-watershot1_80V_1avg_475tgc_1cyc_1.4045MHz_02-Nov-2022 16-18-28.mat';
tic;
disp('Loading watershot data ...');
load([dataset, filesep, filename], 'rcvData', 'settings');
disp(['Completed in ', num2str(toc),  's']);
fprintf('\n');

filename = 'tof_watershot1_02-Nov-2022 16-18-28';
load([dataset, filesep, filename], 'i_tof_water');

% Load mask showing which elements are opposite one another
load([repo_dir, filesep, 'calibration', filesep, 'ideal_element_positions.mat'], 'is_opposite', 'angle_TxNorm_to_RxNorm');

rcvData = removeDCOffset(rcvData);
%%
Ntx = size(rcvData, 1);
Nrx = size(rcvData, 2);
Nt = size(rcvData, 3);

figure;
imagesc( squeeze(rcvData(1,:,301:end)), 0.01*[min(rcvData(:)), max(rcvData(:))]);

% example cross talk figure for TUFFC2022
dt = 1 / settings.acqSamplingFrequency;
t_axis = 0:dt:(Nt - 1) * dt;
tx  = 1;
rxA = 178;
rxB = 183;
figure;
hold on
plot(t_axis*1e6, squeeze(rcvData(tx,rxB,:)), 'r', 'linewidth', 1.5 );
plot(t_axis*1e6, squeeze(rcvData(tx,rxA,:)), 'k', 'linewidth', 1.5 );
ylim([-3000, 3000]);
xlim([116, 136]);
xlabel('Time [\mus]');
ylabel('Voltage [au]');
box on
set(gcf, 'position', [-1204 1337 560 201]);

em_length = 300;

crit_angle   = 90;    % data from receivers within this opening angle from the transmitter will be used
mask         = angle_TxNorm_to_RxNorm <= crit_angle;

% copy the upper triangular to the lower triangular (TOF should be same)
i_tof_water(isnan(i_tof_water)) = 0;
i_tof_water = i_tof_water + i_tof_water';
i_tof_water(i_tof_water == 0) = NaN;
figure;
imagesc(i_tof_water);

Nsig   = 100;
Nxtalk = 80;
Nnoise = 50;
Ngap   = 6;

Nall = Nsig + Nxtalk + Nnoise + (2 * Ngap);
i_sig = Nxtalk + Nnoise + (2 * Ngap) + 1;

Nel = Ntx*Nrx;
extract_data = zeros(Ntx, Nrx, Nall);
noise = zeros(Ntx, Nrx, Nnoise);
sig   = zeros(Ntx, Nrx, Nsig);
xtalk = zeros(Ntx, Nrx, Nxtalk);

i_noise1 = i_sig-Ngap-Nxtalk-Ngap-Nnoise;
i_noise2 = i_sig-Ngap-1-Nxtalk-Ngap;
i_xtalk1 = i_sig-Ngap-Nxtalk;
i_xtalk2 = i_sig-Ngap-1;

for tdx = 1:Ntx
    for rdx = 1:Nrx
        if mask(tdx,rdx)
            data = squeeze( rcvData(tdx, rdx, :) );
            i_fm    = i_tof_water(tdx, rdx);
            i_start = i_fm - (Nxtalk + Nnoise + (2 * Ngap));
            i_end = i_fm + Nsig - 1;
            extract = data(i_start:i_end);
            extract_data(tdx, rdx, :) = extract;
            noise(tdx, rdx, :) = extract(i_noise1:i_noise2);
            sig(tdx, rdx, :)   = extract(i_sig:end);
            xtalk(tdx, rdx, :) = extract(i_xtalk1:i_xtalk2);
        end
    end
end

% This figure shows that the xtalk width is about 50
figure;
plot( squeeze(extract_data(1,150,:))' );
hold on
xline(i_sig);
xline(i_noise1);
xline(i_noise2);
xline(i_xtalk1);
xline(i_xtalk2);

% work out which xtalk signals are valid
thresh = 40;
max_rx = squeeze (max(abs(xtalk), [], 3) );
valid = max_rx > thresh;
figure;
imagesc(valid);
axis image

% find where the xtalk starts for each trace (if at all), set the preceding
% noise to NaN so that xtalk power is not underestimated
for tdx = 1:Ntx
    for rdx = 1:Nrx
        if mask(tdx,rdx)
            data = squeeze( xtalk(tdx, rdx, :) );
            above_t = abs(data) > thresh;
            if sum(above_t) > 1
                i_x = find(above_t, 1, 'first');
                xtalk(tdx, rdx, 1:i_x - 1) = NaN;
            else
                xtalk(tdx, rdx, :) = NaN;
            end
        end
    end
end

figure; imagesc(squeeze(xtalk(1,:,:)));
figure; plot(squeeze(xtalk(1,:,:))');
figure; plot(squeeze(noise(1,:,:))');
figure; plot( squeeze(extract_data(1,:,:))' );
xline(i_sig);
xline(i_noise1);
xline(i_noise2);
xline(i_xtalk1);
xline(i_xtalk2);

Pnoise = zeros(Ntx, Nrx);
Pxtalk = zeros(Ntx, Nrx);
Psig   = zeros(Ntx, Nrx);
SNR = NaN * ones(Ntx, Nrx);
XTALK = NaN * ones(Ntx, Nrx);
% Actually calculate SNR and Xtalk in dB
for tdx = 1:Ntx
    for rdx = 1:Nrx
        if mask(tdx,rdx)
            Pnoise(tdx, rdx) = sum(noise(tdx, rdx, :) .^ 2) / Nnoise;
            xtalk_len = length( ~isnan(xtalk(tdx, rdx, :)) );
            Pxtalk(tdx, rdx) = sum(xtalk(tdx, rdx, :) .^ 2, 'omitnan') / xtalk_len;
            Psig(tdx, rdx) = sum(sig(tdx, rdx, :) .^ 2) / Nsig; 

            SNR(tdx, rdx) = 10 * log10(Psig(tdx,rdx) / Pnoise(tdx, rdx));
            XTALK(tdx, rdx) = 10 * log10(Pxtalk(tdx, rdx) / Psig(tdx, rdx));
        end
    end
end

figure;
histogram(SNR);
title('SNR for all data');

% Extract on-axis only SNR
onaxis_SNR = SNR(:);
is_opposite = is_opposite(:);
onaxis_SNR = onaxis_SNR(is_opposite);
onaxis_SNR = onaxis_SNR(~isnan(onaxis_SNR));

face_colour = [175, 238, 238]/255;

figure;
subplot(3, 2, 1);
h = histogram(onaxis_SNR, 56:0.5:66);
set(h,'facecolor',face_colour);
xlabel('On-axis SNR [dB]');
ylabel('Counts');
xlim([56, 66]);
ylim([0, 50]);

% Vectorise and apply element mask to xtalk
XTALK = XTALK(:);
XTALK = XTALK(mask(:));
XTALK = XTALK(~isinf(XTALK));

figure;
subplot(3, 2, 1);
h = histogram(XTALK, -56:2:-24);
set(h,'facecolor',face_colour);
xlabel('Receive Crosstalk [dB]');
ylabel('Counts');
xlim([-56, -24]);
ylim([0, 2100]);

mean_xtalk = mean(XTALK);
std_xtalk = std(XTALK);
mean_SNR = mean(onaxis_SNR);
std_SNR = std(onaxis_SNR);

disp( ['Mean Xtalk: ', num2str(mean_xtalk), ' dB, std: ', num2str(std_xtalk), ' dB'] );
disp( ['Mean SNR (onaxis): ', num2str(mean_SNR), ' dB, std: ', num2str(std_SNR), ' dB'] );
disp( ['Nmask = ', num2str(sum(mask, 'all', 'omitnan'))] );

example = squeeze(extract_data(1,178,:));
% example(i_noise1:i_noise2) = example(i_noise1:i_noise2) / max(example(i_noise1:i_noise2));
% example(i_xtalk1:i_xtalk2) = example(i_xtalk1:i_xtalk2) / max(example(i_xtalk1:i_xtalk2));
% example(i_sig:end) = example(i_sig:end) / max(example(i_sig:end));
figure;
plot( example );
hold on
xline(i_sig);
xline(i_noise1);
xline(i_noise2);
xline(i_xtalk1);
xline(i_xtalk2);
% on axis SNR histogram and stats
% all xtalk histogram and stats
% time series example image with regions identified
% tune parameters and replot

