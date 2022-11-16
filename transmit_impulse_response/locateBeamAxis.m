function pressure_on_axis = locateBeamAxis(pressure, options)
%LOCATEBEAMAXIS identifies which hydrophone linescan trace corresponds to the beam axis 
%
% DESCRIPTION: 
%      locateBeamAxis estimates where the beam axis is by fitting a
%      quadratic curve to the pressure squared integral of all the linescan
%      traces. The location of the peak of the quadratic (zero gradient) is
%      set to the beam axis position. If the quadratic doesn't have a peak
%      in the measurement window, then the the measurement position with
%      the largest pressure squared integral is assumed to correspond to
%      the beam axis. Optionally, this fitting process can be visualised
%      for the first transmitter in the dataset. Npos: number of linescan
%      measurement positions. Nt: number of time samples. Ntdx: number of
%      transmitters in the dataset.
%
% INPUTS:
%     pressure         - [numeric] input pressure array with size
%                        (Npos, Nt, Ntdx) [Pa]
%
% OPTIONAL INPUTS:
%      Plot            - [boolean] Whether to plot an example of the beam
%                        axis identification process for the first
%                        transmitter in the dataset
%
% OUTPUTS:
%     pressure_on_axis - [numeric] array with size (Nt, Ntdx) containing
%                        the pressure trace at the beam axis for each
%                        transmitter
%
% ABOUT:
%      about       - Morgan Roberts
%      date        - 16th November 2022

arguments
    pressure
    options.Plot
end

% detect size of pressure input
Ntx  = size(pressure, 3);
Nt   = size(pressure, 2);
Npos = size(pressure, 1);

% compute pressure-squared-integral for every position
sum_squared_p = squeeze(sum(pressure.^2, 2));

% Find the beam axis position pos_ba
pos_vec          = 1:Npos;
pressure_on_axis = zeros(Nt, Ntx);
for tdx = 1:Ntx
    % Extract pressure-squared-integral for this transmitter
    ssp = sum_squared_p(:,tdx);
     
    % Fit a quadratic of x-position vs pressure-squared-integral
    [xData, yData] = prepareCurveData(pos_vec, ssp);
    ft             = fittype('poly2');      % Set up quadratic fittype
    fitresult      = fit(xData, yData, ft); % Fit model to data.
    
    % Find position where gradient of quadratic is 0
    pos_ba = round(-fitresult.p2/(2*fitresult.p1)); 
    
    % if the gradient position
    if pos_ba > Npos || pos_ba < 1
        [~,pos_ba] = max(ssp);
        warning(['Not possible to find zero-gradient position for tdx ', num2str(tdx)]);
    end

    % store the pressure on axis
    pressure_on_axis(:,tdx) = pressure(pos_ba, :, tdx);

    % Plot the position location result for element 1 if required
    if options.Plot && tdx == 1
        figure;
        subplot(1, 2, 1);
        plot(pos_vec, ssp);
        hold on;
        plot( fitresult, xData, yData);
        xline(pos_ba, 'k', 'linewidth', 2);
        xlabel('x-position [samples]');
        ylabel('Pressure Squared Integral')

        subplot(1, 2, 2);
        imagesc(squeeze(pressure(:,:,tdx)));
        yline(pos_ba, 'k', 'linewidth', 2);
        ylabel('x-position [samples]');
        xlabel('Time [samples]');
        drawnow
    end

end
end