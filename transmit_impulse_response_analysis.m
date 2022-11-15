% FILENAME:
%     transmit_impulse_response_analysis.m
% 
% DESCRIPTON:
%     Code to analyse the single point hydrophone data measuring the
%     impulse response from several UST transducer elements.
%
% INPUT DATA DIRECTORY:
%     Datasets\open-UST-imaging-and-evaluation\transmit-bandwidth
%
% INPUT DATA FILENAMES:
%     - impulse_response_probe_q.mat      module Q (without matching layers)
%     - impulse_response_probes_AEFG.mat  modules A, E, F and G (with matching layers)
%     - impulse_response_probes_A22uH.mat module A and electrical impedance matching inductors
%     - experimental_parameters.md        details of the experimental setup and data extraction/calibration procedure
%
% FIGURE OUTPUT DIRECTORY:
%     \transmit-bandwidth\figures (cd to local git repo first)
%
% FIGURE OUTPUT FILENAMES:
%       filename  subplot [2,1] (time-vs-voltage and frequency-vs-voltage) showing mean plus measured variation for impulse_response_probe_q
%       filename  subplot [2,1] (time-vs-voltage and frequency-vs-voltage) showing mean plus measured variation for impulse_response_probes_AEFG
%       filename  subplot [2,1] (time-vs-voltage and frequency-vs-voltage) showing mean plus measured variation for impulse_response_probes_A22uH
%       filename  subplot [2,1] (time-vs-voltage and frequency-vs-voltage) showing mean of all three datasets
%       filename  subplot [1,3] histogram of centre-frequency, bandwidth for impulse_response_probe_q
%       filename  subplot [1,3] histogram of centre-frequency, bandwidth for impulse_response_probes_AEFG
%       filename  subplot [1,3] histogram of centre-frequency, bandwidth for impulse_response_probes_A22uH
%       filename  subplot [1,3] histogram of centre-frequency, bandwidth for all three datasets
%
% CALCULATED NUMBERS:
%     - Mean and standard deviation for centre frequeny [Hz]
%     - Mean and standard deviation for amplitude at centre frequency
%     - Mean and standard deviation for bandwidth
%
% ABOUT:
%     - Author:       - Morgan Roberts
%     - Last update:  - 15/11/22

% ---------------------
% Load data
% ---------------------
