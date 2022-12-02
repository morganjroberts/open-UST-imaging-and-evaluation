function [pressure_shift] = alignInTime(pressure, options)
%ALIGNINTIME shifts pressure traces in time to align them
%
% DESCRIPTION:
%      alignInTime uses cross correlation to align pressure traces, and
%      optionally plots the result. Nt: number of time samples. Ntdx:
%      number of transmitters in the dataset. It is assumed that the traces
%      have previously been padded by applyFilterVolume(), and that the
%      padding length is greater than any of the relative delays between
%      traces.
%
% INPUTS:
%      pressure          - [numeric] array with size (Nt, Ntdx) containing
%                          the pressure trace from each transmitter [Pa].
%                          This is expected to be padded previously by
%                          applyFilterVolume().
%
% OPTIONAL INPUTS:
%      Plot              - [boolean] Whether to plot the aligned traces.
%
% OUTPUTS:
%     pressure_shift     - [numeric] array with size (Nt, Ntdx) containing
%                          the aligned pressure traces [Pa]
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

master_p  = pressure(:,1);
pressure_shift = zeros(size(pressure));

for tdx = 1:Ntdx
    % compute cross correlation
    [r,lags] = xcorr(pressure(:,tdx), master_p);
    
    % find lag that gives best alignment
    [~,lagI] = max(r);
    lag      = lags(lagI);
    
    % shift trace
    pressure_shift(:,tdx) = circshift(pressure(:,tdx), -lag);
end

if options.Plot
    figure;
    plot(pressure_shift);
    xlabel('Time [samples]');
    ylabel('Pressure [Pa]');
end

end