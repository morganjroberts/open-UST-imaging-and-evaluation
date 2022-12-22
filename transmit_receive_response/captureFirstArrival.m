function [rcvDataMod, i_fms, i_lms] = captureFirstArrival(rcvData, mask, options)
%CAPTUREFIRSTARRIVAL modifies UST rcvData to contain only the first arrival
%wave packet
%
% DESCRIPTION:
%     captureFirstArrival modifies the receive-data from a UST watershot
%     acquisition. The electromagnetic pickup at ringup, and DC frequency
%     components are both removed. Next, the first-arrival-wave-packet is
%     detected using minimumAIC, and a custom Tukey window is used to set
%     all samples outside of this wave packet to zero, isolating the first
%     arrival. Ntdx: number of transmitters, Nrdx: number of receivers, Nt:
%     number of time samples.
%
% USAGE:
%     rcvData = captureFirstArrival(rcvData, mask)
%
% INPUTS:
%     rcvData       - [numeric] array with size (Ntdx, Nrdx, Nt)
%                     containing watershot acquisition data from a UST
%                     experiment
%     mask          - [boolean] logical array with size (Ntdx, Nrdx)
%                     indicating which transmitter-receiver pairs to
%                     process
%
% OPTIONAL INPUTS:
%     EMlength      - [numeric] integer number of samples required for the
%                     EM pickup to ring down to the noise floor
%     ExtraPlot     - [boolean] whether to plot additional figures showing
%                     the analysis procedure.
%     NoiseLevel    - [numeric] scalar value defining the standard
%                     deviation of the gaussian noise to add to the
%                     signals. Increasing this make the minimumAIC function
%                     less susceptible to pre-first-arrival cross talk.
%     TaperWidth    - [numeric] Integer number of sampels to taper the
%                     signal down to zero, either side of the first/last
%                     motion sample.
%     CaptureWidth  - [numeric] integer number of samples defining the
%                     width of the first-arrival-wave-packet, chosen
%                     visually.
%     PlotIdx       - [numeric] 2 element integer vector [tdx, rdx]
%                     containing the transmitter and receiver index to plot
%
% OUTPUTS:
%     rcvDataMod    - [numeric] array with size (Ntdx, Nrdx, Nt) containing
%                     modified rcvData.
%     i_fms         - [numeric] array with size (Ntdx, Nrdx) containing the
%                     first-motion sample index of the
%                     first-arrival-wave-packet, relative to the start of
%                     the receive data.
%     i_lms         - [numeric] array with size (Ntdx, Nrdx) containing the
%                     last-motion sample index of the
%                     first-arrival-wave-packet, relative to the start of
%                     the receive data.
%
% DEPENDENCIES:
%     minimumAIC / pickFirstMotionAIC
%
% ABOUT:
%     author        - Morgan Roberts
%     date          - 21/11/2022

arguments
    rcvData
    mask
    options.EMlength      = 300;
    options.ExtraPlot     = true;
    options.NoiseLevel    = 1e-2;
    options.TaperWidth    = 6;
    options.CaptureWidth  = 150;
    options.PlotIdx       = [1,144];
end

% set random seed to default for repeatable results
rng('default');

% calculate number of transmitters/receivers
Ntdx = size(rcvData, 1);
Nrdx = size(rcvData, 2);
Nt   = size(rcvData, 3);

% remove DC offset
rcvData = removeDCOffset(rcvData);

% replace the electromagnetic interference region at the start with zeros
rcvData(:,:,1:options.EMlength) = 0;

% Initialise storage for outputs
i_fms      = NaN * ones(Ntdx, Nrdx);
i_lms      = NaN * ones(Ntdx, Nrdx);
rcvDataMod = NaN * ones(size(rcvData));

% Progress tracing
counter    = 0;
prog_round = -1;

tic;
disp('Capturing first arrivals ...');
for tdx = 1:Ntdx
    progress = floor(1e2 * (counter / Ntdx) / 10) * 10;
    if progress > prog_round
        disp(['Progress: ', num2str(progress), ' % after ', num2str(toc), ' s'])
        prog_round = progress;
    end
    counter = counter + 1;

    for rdx = 1:Nrdx
        if mask(tdx, rdx)
            % Extract signal
            trace = squeeze(rcvData(tdx, rdx, :));

            % Reset random seed, add gaussian noise to help the AIC pick
            rng(1);
            noise   = options.NoiseLevel * randn(size(trace));
            trace_n = noise + (trace / max(trace));

            % Find index of peak value in initial wave packet (to stop AIC picking after the main peak)
            [~, i_max] = max(trace_n(:));

            % Detect first-motion-sample using AIC and store
            i_fm            = minimumAIC(trace_n(1:i_max));
            i_fms(tdx, rdx) = i_fm;

            if options.ExtraPlot && tdx == 1 && rdx == 206
                minimumAIC(trace_n(1:i_max), Plot=true);
            end

            % Build a custom tukey window for the signal ring up
            win   = getWin((2 * options.TaperWidth), 'Tukey', 'Param', 1);
            win_L = [zeros(i_fm - 1 - options.TaperWidth, 1); win(1:options.TaperWidth)];
            
            % Apply the window to all samples before the first-motion
            trace(1:i_fm-1) = trace(1:i_fm-1) .* win_L;
%             rcvDataMod(tdx, rdx, 1:i_fm-1) = squeeze(rcvData(tdx, rdx, 1:i_fm-1)) .* win_L;

            % Find last-motion-sample based on a user-defined capture width
            i_lm            = i_fm + options.CaptureWidth;
            i_lms(tdx, rdx) = i_lm;

            % Build a custom tukey window for the signal ring down
            win_R = [(options.TaperWidth+1:end); zeros(Nt - i_lm - options.TaperWidth, 1)];
            
            % Apply the window to all samples after the last-motion
            trace(i_lm + 1:end)                = trace(i_lm + 1:end) .* win_R;
%             rcvDataMod(tdx, rdx, i_lm + 1:end) = squeeze(rcvData(tdx, rdx, i_lm + 1:end)) .* win_R;

            % Save the modified rcvData
            rcvDataMod(tdx, rdx, :) = trace;

            if options.ExtraPlot && tdx == options.PlotIdx(1) && rdx == options.PlotIdx(2)
                figure;
                subplot(2, 1, 1);
                hold on
                plot(squeeze(rcvData(tdx, rdx, :)), 'k')
                xline(i_fm, 'k--')
                xline(i_lm, 'k--')
                xlim([i_fm-100, i_lm+100]);
                xlabel('Time [samples]');
                ylabel('Voltage [au]');
                title('Raw Data');

                subplot(2, 1, 2);
                hold on
                plot(trace, 'k')
                xline(i_fm, 'k--')
                xline(i_lm, 'k--')
                xlim([i_fm-100, i_lm+100]);
                xlabel('Time [samples]');
                ylabel('Voltage [au]');        
                title(['Processed Trace, tx = ', num2str(options.PlotIdx(1)), ', rx = ', num2str(options.PlotIdx(2))]);

                drawnow    
            end  
        end
    end

    if options.ExtraPlot && tdx == options.PlotIdx(1)
        figure;
        hold on;
        for rdx = 1:Nrdx
            if mask(1, rdx)
                rx_data = squeeze(rcvData(1,rdx,:));
                rx_data = rx_data / max(abs(rx_data));
                plot(rx_data + rdx, 1:Nt, 'k');
            end 
        end
        xlabel('Receivers');
        ylabel('Time [samples]');
        h1 = plot(i_fms(1,:), 'r');
        h1.MarkerEdgeColor= 'r';
        h1.Marker = '.';
        h2 = plot(i_lms(1,:), 'r');
        h2.MarkerEdgeColor= 'r';
        h2.Marker = '.';    
        ylim([min(i_fms, [], 'all', 'omitnan'), max(i_lms, [], 'all', 'omitnan')]);
        title(['First arrival capture, tx = ', num2str(options.PlotIdx(1))]);
        drawnow

    end 

end
disp(['Completed in ', num2str(toc),  's']);
fprintf('\n');


end