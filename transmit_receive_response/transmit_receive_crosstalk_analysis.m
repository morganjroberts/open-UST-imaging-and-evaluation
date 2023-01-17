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
load([repo_dir, filesep, 'calibration', filesep, 'ideal_element_positions.mat'], 'is_opposite');

rcvData = removeDCOffset(rcvData);
%%
Ntx = size(rcvData, 1);
Nrx = size(rcvData, 2);
Nt = size(rcvData, 3);



% replace EM pickup with noise
% noise_sample = squeeze( rcvData(1, 178, 301:600) );
% noise_sample = permute(noise_sample, [2, 3, 1]);
% noise_sample = repmat(noise_sample, [Ntx, Nrx]);
% rcvData(:,:,1:300) = noise_sample;

em_length = 300;

element_mask = ~isnan(i_tof_water);

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
% SNR = zeros(Ns;

i_noise1 = i_sig-Ngap-Nxtalk-Ngap-Nnoise;
i_noise2 = i_sig-Ngap-1-Nxtalk-Ngap;
i_xtalk1 = i_sig-Ngap-Nxtalk;
i_xtalk2 = i_sig-Ngap-1;

for tdx = 1:Ntx
    for rdx = 1:Nrx
        if element_mask(tdx,rdx)
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
        if element_mask(tdx,rdx)
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
        if element_mask(tdx,rdx)
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

figure;
histogram(onaxis_SNR);
xlabel('SNR [dB]');
ylabel('Counts');

% Vectorise and apply element mask to xtalk
XTALK = XTALK(:);
XTALK = XTALK(element_mask(:));
XTALK = XTALK(~isinf(XTALK));
figure;
histogram(XTALK)
xlabel('Receive Crosstalk [dB]');
ylabel('Counts');

mean_xtalk = mean(XTALK);
std_xtalk = std(XTALK);
mean_SNR = mean(onaxis_SNR);
std_SNR = std(onaxis_SNR);

disp( ['Mean Xtalk: ', num2str(mean_xtalk), ' dB, std: ', num2str(std_xtalk), ' dB'] );
disp( ['Mean SNR (onaxis): ', num2str(mean_SNR), ' dB, std: ', num2str(std_SNR), ' dB'] );

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

