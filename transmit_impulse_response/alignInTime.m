function [pressure_shift_pad, pressure_shift] = alignInTime(pressure, options)
%ALIGNINTIME shifts pressure traces in time to align them
%
% DESCRIPTION:
%      alignInTime uses cross correlation to align pressure traces, and
%      optionally plots the result. Nt: number of time samples. Ntdx:
%      number of transmitters in the dataset. It is assumed that the traces
%      have previously been padded by filterTraces(), and that the padding
%      length is greater than any of the relative delays between traces.
%      Two versions of the output array are provided: the first with the
%      padding intact, to be fed into computeSpect(), and the second with
%      the padding removed, to be fed into plotMeanAndVariation().
%
% INPUTS:
%      pressure          - [numeric] array with size (Nt, Ntdx) containing
%                          the pressure trace from each transmitter [Pa].
%                          This is expected to be padded previously by
%                          filterTraces().
%
% OPTIONAL INPUTS:
%      Plot              - [boolean] Whether to plot the aligned traces.
%
% OUTPUTS:
%     pressure_shift_pad - [numeric] array with size (Nt, Ntdx) containing
%                          the aligned pressure traces, with the padding
%                          intact [Pa]
%     pressure_shift     - [numeric] array with size (Nt/3, Ntdx) containing
%                          the aligned pressure traces, with the padding
%                          removed [Pa]
%
% ABOUT:
%      author            - Morgan Roberts
%      date              - 16th November 2022


arguments
    pressure
    options.Plot
end

% detect size of pressure input
Ntdx = size(pressure, 2);
Nt   = size(pressure, 1);

master_p  = pressure(:,1);
pressure_shift_pad = zeros(size(pressure));

for tdx = 1:Ntdx
    % compute cross correlation
    [r,lags] = xcorr(pressure(:,tdx), master_p);
    
    % find lag that gives best alignment
    [~,lagI] = max(r);
    lag      = lags(lagI);
    
    % shift trace
    pressure_shift_pad(:,tdx) = circshift(pressure(:,tdx), -lag);
end

if options.Plot
    figure;
    plot(pressure_shift_pad);
    xlabel('Time [samples]');
    ylabel('Pressure [Pa]');
end

% Remove the padding for the next function
pressure_shift = pressure_shift_pad( (Nt / 3) + 1 : end - (Nt / 3) , :);

end