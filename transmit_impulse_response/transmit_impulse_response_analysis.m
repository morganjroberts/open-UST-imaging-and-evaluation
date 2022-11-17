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
%     transmit_impulse_response_probe_q.mat        module Q (without acoustic matching layers)
%     transmit_impulse_response_probes_AEFG.mat    modules A, E, F and G (with acoustic matching layers)
%     transmit_impulse_response_probes_A22uH.mat   module A driven thorugh electrical impedance matching inductors
%
% EXPERIMENTAL PARAMETERS:
%     \transmit_impulse_response\transmit_impulse_response_experimental_parameters.md (cd to local git repo first)
%
% FIGURE OUTPUT DIRECTORY:
%     \transmit_impulse_response\figures (cd to local git repo first)
%
% CALCULATED NUMBERS:
%     - Mean and standard deviation for centre frequeny [Hz]
%     - Mean and standard deviation for amplitude at centre frequency
%     - Mean and standard deviation for bandwidth
%
% ABOUT:
%     author:      - Morgan Roberts
%     last update: - 15/11/22

close all
clearvars
warning('off')

% ---------------------
% Load data
% ---------------------
% change this to match the correct mapped drive location
input_dir = 'Z:\open-UST-imaging-and-evaluation\transmit_impulse_response\';

% filtering cutoff frequency
cutoff_f = 5e6;

% load the files and process the data
% input_filename = 'transmit_impulse_response_probe_q';
% [t_Q, mean_p_Q, f_Q, mean_as_Q] = processImpulseResponse(input_dir, input_filename, ...
%     cutoff_f, ExtraPlot=1, spectMode='db', dbThresh=-10, traceXlim='time_axis', traceYlim=[-Inf, Inf], spectXlim=[0, 5], spectYlim=[-40,2]);

% input_filename = 'transmit_impulse_response_probe_A22uH';
% [t_A22uH, mean_p_A22uH, f_A22uH, mean_as_A22uH] = processImpulseResponse(input_dir, input_filename, ...
%     cutoff_f, ExtraPlot=0, spectMode='db', dbThresh=-10, traceXlim='time_axis', traceYlim=[-Inf, Inf], spectXlim=[0, 5], spectYlim=[-40,2]);
%  
input_filename = 'transmit_impulse_response_probes_AEFG';
[t_AEFG, mean_p_AEFG, f_AEFG, mean_as_AEFG] = processImpulseResponse(input_dir, input_filename, ...
    cutoff_f, ExtraPlot=0, spectMode='db', dbThresh=-10, traceXlim='time_axis', traceYlim=[-Inf, Inf], spectXlim=[0.4, 5], spectYlim=[-40,2]);

%% Compare the three datasets

figure;
subplot(2, 1, 1);
hold on;
plot(t_Q*1e6, mean_p_Q*1e-3)
plot(t_AEFG*1e6, mean_p_AEFG*1e-3)
plot(t_A22uH*1e6, mean_p_A22uH*1e-3)
xlabel('Time [us]');
ylabel('Pressure [kPa]');

subplot(2, 1, 2);
hold on;
plot(f_Q*1e-6, mean_as_Q)
plot(f_AEFG*1e-6, mean_as_AEFG)
plot(f_A22uH*1e-6, mean_as_A22uH)
xlim([0, 5])
xlabel('Frequency [MHz]');
ylabel('Amplitude [Pa]');
legend({'No Matching (N=16)', 'Acoustic Matching Only (N=64)', 'Acoustic and Electrical Matching (N=16)'});

% % Compare the matched and unmatched probes
% figure;
% f_qu       = linspace(0, 5e6, 1000);
% matched    = interp1(f_AEFG, mean_as_AEFG, f_qu);
% un_matched = interp1(f_Q, mean_as_Q, f_qu);

