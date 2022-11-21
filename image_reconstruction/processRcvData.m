function rcvData = processRcvData(rcvData, options)

% EMlength - length of electromagnetic interference at start of signal to
% be replaced by zeros
% Symmetric - boolean stating whether to average the reciprocal tx/rx
% waveforms to make the rcvData syymetric

arguments
    rcvData (:,:,:) {mustBeNumeric, mustBeReal};
    options.EMlength (1,1) {mustBeInteger, mustBeReal} = 300;
    options.Symmetric (1,1) logical = true;
end

Ntx = size(rcvData, 1);
Nrx = size(rcvData, 2);

% remove DC offset
rcvData = removeDCOffset(rcvData);

% replace the electromagnetic interference region at the start with zeros
rcvData(:,:,1:options.EMlength) = 0;

% every tx/rx pair has two signals : tx==>rx and rx==>tx. This part of the
% code takes the mean of both signals.
% this makes rcvData symmetric for all time points

if options.Symmetric
    for tdx = 1:Ntx
        for rdx = 1:Nrx
            if tdx < rdx
                % extract similar reciprocal traces
                txrx_trace = squeeze(rcvData(tdx, rdx, :));
                rxtx_trace = squeeze(rcvData(rdx, tdx, :));    
                
                % calculate mean
                mean_trace = mean([txrx_trace, rxtx_trace], 2);
                
                % reassign identical traces to tx/rx pair
                rcvData(tdx,rdx,:) = mean_trace;
                rcvData(rdx,tdx,:) = mean_trace;
            end
        end
    end
end