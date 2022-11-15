% FILENAME:
%     transmit_impulse_response_analysis.m
% 
% DESCRIPTON:
%     Code to analyse the single point hydrophone data measuring the
%     impulse response from several UST transducer elements. First the data
%     is filtered to remove the high frequency and DC content. Next, the
%     x-position corresponding to the beam axis is located for every
%     element. All other non-beam-axis traces are then discarded. Next, the
%     waveforms are windowed, pre/post padded, and aligned in time. The
%     amplitude spectra are then calculated using an FFT. Finally,
%     everything is plotted, and relevant quantities are calculated.
%
% INPUT DATA DIRECTORY:
%     Datasets\open-UST-imaging-and-evaluation\transmit_impulse_response
%
% INPUT DATA FILENAMES:
%     transmit_impulse_response_probe_q.mat                module Q (without acoustic matching layers)
%     transmit_impulse_response_probes_AEFG.mat            modules A, E, F and G (with acoustic matching layers)
%     transmit_impulse_response_probes_A22uH.mat           module A driven thorugh electrical impedance matching inductors
%
% EXPERIMENTAL PARAMETERS:
%     transmit_impulse_response_experimental_parameters.md 
%
% FIGURE OUTPUT DIRECTORY:
%     \transmit-impulse_response\figures (cd to local git repo first)
%
% FIGURE OUTPUT FILENAMES:
%     filename  subplot [2,1] (time-vs-voltage and frequency-vs-voltage) showing mean plus measured variation for impulse_response_probe_q
%     filename  subplot [2,1] (time-vs-voltage and frequency-vs-voltage) showing mean plus measured variation for impulse_response_probes_AEFG
%     filename  subplot [2,1] (time-vs-voltage and frequency-vs-voltage) showing mean plus measured variation for impulse_response_probes_A22uH
%     filename  subplot [2,1] (time-vs-voltage and frequency-vs-voltage) showing mean of all three datasets
%     filename  subplot [1,3] histogram of centre-frequency, bandwidth for impulse_response_probe_q
%     filename  subplot [1,3] histogram of centre-frequency, bandwidth for impulse_response_probes_AEFG
%     filename  subplot [1,3] histogram of centre-frequency, bandwidth for impulse_response_probes_A22uH
%     filename  subplot [1,3] histogram of centre-frequency, bandwidth for all three datasets
%
% CALCULATED NUMBERS:
%     - Mean and standard deviation for centre frequeny [Hz]
%     - Mean and standard deviation for amplitude at centre frequency
%     - Mean and standard deviation for bandwidth
%
% ABOUT:
%     author:       - Morgan Roberts
%     last update:  - 15/11/22

close all
clearvars

% ---------------------
% Load data
% ---------------------
% change this to match the correct mapped drive location
input_dir = 'Z:\open-UST-imaging-and-evaluation\transmit-bandwidth\';


function voltage_filt = filterTraces(voltage, dt, filter_param)



end

