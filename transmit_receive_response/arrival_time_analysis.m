% FILENAME:
%     .m
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

% copy the upper triangular to the lower triangular (TOF should be same)
i_tof_water(isnan(i_tof_water)) = 0;
i_tof_water = i_tof_water + i_tof_water';
i_tof_water(i_tof_water == 0) = NaN;

dt = 1 / settings.acqSamplingFrequency;

tof = i_tof_water(:) * dt;
tof = tof(is_opposite(:));

mean_tof = mean(tof);
std_tof = std(tof);
std_d = std_tof * 1489;
std_N = std_tof / (1 / 1.21e6);
disp( ['Mean TOF: ', num2str(1e6*mean_tof), ' us, std: ', num2str(1e6*std_tof), ' us'] );
disp( ['Std position error: ', num2str(std_d*1e3), 'mm'] );
disp( ['Std cycles error: ', num2str(std_N), 'cyc'] );

figure;
histogram(tof);