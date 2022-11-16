function [freqs, as] = computeSpect(pressure, dt)
%COMPUTESPECT applies a Tukey window and calculates amplitude spectrum
%
% DESCRIPTION:
%     computeSpect takes a pressure trace dataset that has already been
%     pre/post zero-padded, applies a sharp tukey window, and performs a
%     FFT to calculate the amplitude spectrum. Nt: number of time samples.
%     Ntdx: number of transmitters in the dataset.
%
% INPUTS:
%      pressure      - [numeric] pressure array with size (Nt, Ntdx) [Pa]
%      dt            - [numeric] scalar time step [s]
%
% OUTPUTS:
%     freqs          - [numeric] frequency vector with Nf frequencies [Hz]
%     as             - [numeric] array with size (Nf, Ntdx) containing the
%                      amplitude spectrum for each transmitter [Pa]
%
% ABOUT:
%      author        - Morgan Roberts
%       date         - 16th November 2022

% detect size of pressure input
Ntdx = size(pressure, 2);
Nt   = size(pressure, 1);

% calculate sampling frequency [Hz]
Fs = 1 / dt;

% apply a sharp Tukey window
win          = getWin(Nt, 'Tukey', 'Param', 0.05);
win          = repmat(win, [1, Ntdx]);
pressure_win = win .* pressure;

% Compute the amplitude spectrum
[freqs, as] = spect(pressure_win, Fs, 'FFTLength', Nt*6, 'Dim', 1);
end