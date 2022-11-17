% INPUT DATA DIRECTORY:
%     Datasets\open-UST-imaging-and-evaluation\UST_dataset1
%
% INPUT DATA FILENAMES:
%     <input-data-dir>\UST_acquisition_UST-watershot1_80V_1avg_475tgc_1cyc_1.4045MHz_02-Nov-2022 16-18-28.mat  UST watershot dataset

close all
clearvars

data_dir = 'Z:\open-UST-imaging-and-evaluation\UST_dataset1\';
filename = 'UST_acquisition_UST-watershot1_80V_1avg_475tgc_1cyc_1.4045MHz_02-Nov-2022 16-18-28.mat';
load([data_dir, filename], 'rcvData', 'settings', 'temperature', 'time_axis', 'wvfmF');

