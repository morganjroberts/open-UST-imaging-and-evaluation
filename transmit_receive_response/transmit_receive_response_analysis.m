% FILENAME:
%     transmit_receive_response_analysis.m
% 
% DESCRIPTON:
%      Code to analyse a UST watershot dataset from the open-UST transducer
%      array. For each transmit event, a single transmitter generates an
%      ultrasound wave and all other elements act as recievers. Every
%      transmit event uses a different transmitting element. Each trace in
%      the watershot dataset is the receive data from a unique
%      transmitter/receiver pair. The transmit-receive response is analysed
%      to assess the effect of interelement variation on time-of-arrival,
%      and amplitude spectrum content, and also to assess the effect of
%      directionality on these metrics.
%
% INPUT DATA FILENAMES:
%     <data-dir>\UST_dataset1\UST_acquisition_UST-watershot1_80V_1avg_475tgc_1cyc_1.4045MHz_02-Nov-2022 16-18-28.mat  UST watershot dataset
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
[repo_dir,data_dir] = getRepoDataPath;
dataset  = [data_dir, filesep, 'UST_dataset1'];
filename = 'UST_acquisition_UST-watershot1_80V_1avg_475tgc_1cyc_1.4045MHz_02-Nov-2022 16-18-28.mat';
tic;
disp('Loading watershot data ...');
load([dataset, filesep, filename], 'rcvData', 'settings', 'temperature', 'time_axis', 'wvfmF');
disp(['Completed in ', num2str(toc),  's']);
fprintf('\n');

% Load predicted element positions
load([repo_dir, filesep, 'calibration', filesep, 'ideal_element_positions.mat'], 'distance_matrix', 'angle_to_norm_matrix');

% Experimental parameters
temp = 21.7; % measurement temp
Fs   = settings.acqSamplingFrequency;

% Create a mask
opening_angle = 1;    % data from receivers within this opening angle from the transmitter will be used
mask          = angle_to_norm_matrix <= opening_angle / 2;

% Capture the first arrivals of the pulses
fa_data = captureFirstArrival(rcvData, mask, ExtraPlot=true, NoiseLevel=2.5e-2, CaptureWidth=150, WinWidth=6);

% Compute the amplitude spectra
Ntdx    = size(fa_data, 1);
Nrdx    = size(fa_data, 2);
Nt      = size(fa_data, 3);
fa_data = reshape(fa_data, [Ntdx*Nrdx,Nt]);
fa_data = fa_data(mask(:),:);
[f, as] = spect(fa_data, Fs, 'Dim', 2);
figure;
plot(fa_data')



% group pulses into angle-bins

% for the on-axis bin, look at time shifting - if it's not good enough,
% align in time e.g. xcorr, work out distribution of lags

% plot on-axis time and amplitude spect

% plot angle vs amplitude spect
