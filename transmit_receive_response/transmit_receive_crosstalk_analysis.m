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

sig_width   = 130;
xtalk_width = 50;

pre_arrival = zeros(size(rcvData));
SNR = zeros(Ntx, Nrx);

for tdx = 1:Ntx
    for rdx = 1:Nrx
        if element_mask(tdx,rdx)
            data = squeeze( rcvData(tdx, rdx, :) );
            pre_data = data(em_length:i_tof_water(tdx, rdx));
            N = length(pre_data);
            pre_arrival(tdx, rdx, end-N+1:end) = pre_data;
        end
    end
end

% This figure shows that the xtalk width is about 50
figure;
plot( squeeze(pre_arrival(1,:,:))' );

% An example trace
figure;
hold on
plot( squeeze(rcvData(1, 128, :)) )
xline(i_tof_water(1, 128));

% Need to detect which channels have cross talk by removing DC and using
% threshold
thresh_val = 50;

for tdx = 1:Ntx
    for rdx = 128%:Nrx
        if element_mask(tdx,rdx)
            data   = squeeze( rcvData(tdx, rdx, :) );
            figure;
            plot(data);
            keyboard
        end
    end
end
keyboard

% actually compute SNR
for tdx = 1:Ntx
    for rdx = 1:Nrx
        if element_mask(tdx,rdx)
            data   = squeeze( rcvData(tdx, rdx, :) );
            thresh = mean(data) + thresh_val;
            if max(abs(pre_arrival(tdx, rdx, :))) > thresh
%                 disp('here');
                i_c    = i_tof_water(tdx, rdx);
                xtalk  = data( i_c-xtalk_width : i_c );
                sig    = data( i_c : i_c+sig_width );
                Pxtalk = sum( xtalk .^ 2 ) / length(xtalk);
                Psig   = sum( sig .^ 2 ) / length(sig);
                SNR(tdx, rdx) = 10 * log10( Pxtalk / Psig );
            end
        end
    end
end

% exract the non-NaN values
SNR = SNR(:);
mask_vec = (SNR ~= 0);
SNR = SNR(mask_vec);

figure;
histogram(SNR);

