function g_data = groupSimilarSignals(data, mask, options)
%GROUPSIMILARSIGNALS extracts, aligns and reshapes pre-modified UST rcvData according to mask
%
% DESCRIPTION:
%     groupSimilarSignals takes pre-modified UST rcvData (e.g. from
%     captureFirstArrival) and extracts the signals corresponding to tx/rx
%     pairs defined by a selection mask. The selection mask can be used to
%     select similar signals, e.g. those rays with a similar
%     emission/incidence angle. Optionally the signals can be aligned using
%     cross correlation beofre being returned. To prevent cycle skipping errors,
%     a maximum delay can be specified for the alignment process. The
%     rcvData has Ntdx transmitters, Nrdx recievers, and Nt time points.
%
% USAGE:
%      g_data = groupSimilarSignals(data, mask)
%
% INPUTS:
%     data          - [numeric] array with size (Ntdx, Nrdx, Nt)
%                     containing pre-modified UST rcvData
%     mask          - [boolean] logical array with size (Ntdx, Nrdx)
%                     indicating which transmitter-receiver pairs to
%                     process
%
% OPTIONAL INPUTS:
%     Align         - [boolean] Whether to align the signals using cross
%                     correlation before returning them
%     RestrictXcorr - [boolean] maximum amount of delay
%                     permitted between the first-motion sample of the
%                     signals, after alignment
%     ExtraPlot     - [boolean] whether to plot additional figures showing
%                     the analysis procedure.
%
% OUTPUTS:
%     g_data        - [numeric] grouped data with size (Nel, Nt) where Nel
%                     is the number of valid tx/rx pairs in the mask.
%                     Returns a scalar value of zero if the mask is empty.
%
% ABOUT:
%     author        - Morgan Roberts
%     date          - 23/11/2022

arguments
    data
    mask
    options.Align         = true;
    options.RestrictXcorr = 6;
    options.ExtraPlot     = true
end

% Return a scalar value of zero as an error flag if the mask is empty
if sum(mask) < 1
    g_data = 0;
    return
end

% Compute sizes of array
Ntdx = size(data, 1);
Nrdx = size(data, 2);
Nt   = size(data, 3);

% Reshape the data since original tx/rx info no longer important
data  = reshape(data, [Ntdx * Nrdx, Nt]);
mask  = reshape(mask, [Ntdx * Nrdx, 1]);

% Extract the grouped data
g_data = squeeze(data(mask, :));

% Align the grouped data in time if required
if options.Align
    g_data_aligned = zeros(size(g_data));
    
    N = size(g_data, 1);
    
    master = g_data(1, :);
    tic;
    disp('Grouping similar signals ...');
    for idx = 1:N
        % Make copy of the current trace
        trace = squeeze(g_data(idx,:));
    
        % Compute cross correlation
        [r, lags] = xcorr(master, trace);
    
        % Modify the cross correltaion result to exclude lag positions outside
        % of expected range (estimated from the delay between first-motion of each)
        del     = find(master, 1, 'first') - find(trace, 1, 'first');
        lag_min = del - options.RestrictXcorr;
        lag_max = del + options.RestrictXcorr;    
        r(1:find(lags == lag_min) - 1)     = 0;
        r(find(lags == lag_max) + 1 : end) = 0;
    
        % Find the best alignment and store
        [~,I] = max(r);
        lag   = lags(I);
        g_data_aligned(idx, :) = circshift(trace, lag);
        
    end
    disp(['Completed in ', num2str(toc),  's']);
    fprintf('\n');
end

% Plot the data if required
if options.ExtraPlot
    i_start = find(sum(g_data, 1), 1, 'first') - 10;
    i_end   = find(sum(g_data, 1), 1, 'last') + 10;

    figure;

    if options.Align
        subplot(2, 1, 2);
        plot(g_data_aligned')
        title('Aligned')
        xlim([i_start, i_end]);
        xlabel('Time [samples]');
        ylabel('Voltage [au]');

        subplot(2, 1, 1);
    end

    plot(g_data')
    title('Original')
    xlim([i_start, i_end]);
    xlabel('Time [samples]');
    ylabel('Voltage [au]');
    
    drawnow
end

if options.Align
    g_data = g_data_aligned;
end

end