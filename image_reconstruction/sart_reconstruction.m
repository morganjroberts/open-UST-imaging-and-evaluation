% FILENAME:
%     sart_reconstruction.m
% 
% DESCRIPTON:
%     Reconstruct ultrasound tomography (UST) data acquired from 2D ring
%     array. The data is the difference in time-of-flight between a phantom
%     dataset and a waterhost dataset. An estimate of the sound speed
%     distribution in the phantom is reconstructed using straight ray
%     simultaneous-algebraic-reconstruction-technique (SART).
%
% INPUT DATA DIRECTORY:
%     Datasets\open-UST-imaging-and-evaluation\UST_dataset1
%
% INPUT DATA FILENAMES:
%     <input-data-dir>\delta_tof_phantom.mat               time-of-flight picked data from a UST experiment
%     <repo-dir>\calibration\ideal_element_positions.mat   ideal cartesian element positions (X, Y)
%     
% FIGURE OUTPUT DIRECTORY:
%     <repo-dir>\image_reconstruction\figures
%
% DEPENDENCIES:
%     ust-sart https://github.com/ucl-bug/ust-sart
%     k-Wave   https://github.com/ucl-bug/k-wave
%
% ABOUT:
%     author:      - Morgan Roberts
%     last update: - 17/11/22

close all
clearvars

% make sure the data_dir is correct. Current folder should be the repo_dir
data_dir = 'Z:\open-UST-imaging-and-evaluation\UST_dataset1\';
load([data_dir, 'delta_tof_phantom.mat'], 'delta_tof');
load('..\calibration\ideal_element_positions.mat', 'element_positions');

% apply a 3x3 2D median filter to the input data to remove noise
filt_delta_tof = medfilt2(delta_tof);

% Initialise sart object
temperature = 20;   % water temperature [degC]
sart        = SartExperiment(element_positions, filt_delta_tof, temperature);

% Plot setup and data
sart.plotSetup;

% Perform reconstruction
ups        = [1 * ones(1, 22), ...
              2 * ones(1, 17), ...
              4 * ones(1, 15)]; % upsampling factors for each iteration
Nit        = length(ups);  % number of iterations
dx0        = 4e-3;         % step size for iteration 1 [m]
init_c_val = sart.c_water; % sound speed value for homogeneous initial estimate [m/s]
hamming    = 0;            % boolean controlling whether hamming window is used
recon_d    = 0.130;        % diameter of reconstruction circle [m]
sart.reconstructSart(init_c_val, recon_d, Nit, dx0, ups, hamming=hamming);

% Plot and save the reconstruction result
sart.saveReconResult;
sart.plotReconResult;    

