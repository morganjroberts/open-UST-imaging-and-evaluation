% FILENAME:
%     extract_first_arrival_wave_packet.m
% 
% DESCRIPTON:
%      Code to post-process a UST watershot dataset from the open-UST
%      transducer array. The first-arrival-wave-packet is isolated, and all
%      EM pickup and structural reflections are removed. The DC offset is
%      also removed. This extraction process is in a separate script
%      because of its long computation time.
%
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

% Load predicted element positions
load([repo_dir, filesep, 'calibration', filesep, 'ideal_element_positions.mat'], 'element_positions', 'distance_matrix', 'angle_TxNorm_to_RxNorm');

% Create a mask
crit_angle    = 90;    % data from receivers within this opening angle from the transmitter will be used
mask          = angle_TxNorm_to_RxNorm <= crit_angle;
plotActiveElements(element_positions, mask, Tdx=10, Loop=false)

% Capture the first arrivals of the pulses
fa_data = captureFirstArrival(rcvData, mask, ExtraPlot=true, NoiseLevel=2.5e-2, CaptureWidth=130, TaperWidth=8, PlotIdx=[1,70]);

output_filename = 'First_arrival_capture_UST_watershot_2-11-22 16-18-28';
save([dataset, filesep, output_filename], 'fa_data', 'mask', 'settings');
