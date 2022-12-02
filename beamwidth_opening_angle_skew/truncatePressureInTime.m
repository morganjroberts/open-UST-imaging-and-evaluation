function [pressure_win, time_axis_cut, Nt] = truncatePressureInTime(pressure, time_axis, t_end, options)
% TRUNCATEPRESSUREINTIME Truncate a pressure measurement plane in time
%
% DESCRIPTION:
%     truncatePressureInTime is used to truncate the end of a pressure
%     measurement plane, for example to remove a hydrophone-tip reflection.
%     The cut is smoothed with a right-sided Tukey window, and optionally
%     the traces can be plotted before and after truncation to help choose
%     the most suitable cutting location.
%
% USAGE:
%     [pressure_win, time_axis_cut, Nt] = truncatePressureInTime(pressure,
%     time_axis, t_end)
%
% INPUTS:
%     pressure      - [numeric] pressure array with size (Ny, Nx, Nt) [Pa]
%     time_axis     - [numeric] time axis vector with size (1, Nt) [s]
%     t_end         - [numeric] the new end time of the pressure data [us]
%
% OPTIONAL INPUTS:
%     Plot          - [boolean] whether to plot the traces before and after
%                     truncation. 
%     x             - [integer] index for which x-position to plot (column)
%     y             - [integer] index for which y-position to plot (row)
%     Clip          - [numeric] scalar value defining the clip factor for
%                     the colour scale, relative to the maximum value in
%                     the measurement plane
%     TukeyParam    - [numeric] scalar value between 0 and 1 defining the
%                     taper width of the tukey window used to smooth the
%                     truncation cut.
%
% OUTPUTS:
%     pressure_win  - [numeric] truncated pressure array [Pa]
%     time_axis_cut - [numeric] truncated time axis vector [s]
%     Nt            - [integer] new length of the time axis [pts]
%
% ABOUT:
%     author        - Morgan Roberts
%     date:         - 2/12/22

arguments
    pressure
    time_axis
    t_end
    options.Plot       = true;
    options.x          = round(size(pressure, 2) / 2);
    options.y          = round(size(pressure, 1) / 2);
    options.Clip       = 0.02;
    options.TukeyParam = 0.05;
end

% Shorten variable name
ix = options.x;
iy = options.y;

tic;
disp('Removing Hydrophone Tip Reflection...');

% Compute length of spatial dimensions
Nx = size(pressure, 2);
Ny = size(pressure, 1);

% Find exact discretised cutting position
[t_cut, Nt] = findClosest(time_axis, t_end * 1e-6);

% Trim the pressure data and update other parameters
pressure_cut  = pressure(:, :, 1:Nt);
time_axis_cut = time_axis(1:Nt);

% Create a right-sided sided window to smooth the cut, replicate in 3D
win                            = getWin(Nt, 'Tukey', 'Param', options.TukeyParam);
win( 1:find(win == 1, 1) - 1 ) = 1;
win                            = permute(win, [2, 3, 1]);
win                            = repmat(win, [Ny, Nx, 1]);

% Apply the window
pressure_win = pressure_cut .* win;

% Plot the data before and after truncation
if options.Plot
    % extract a slice and trace at the specified y and x positions
    slice_y     = squeeze(pressure(iy, :, :));
    trace       = squeeze(pressure(iy, ix, :));
    slice_y_cut = squeeze(pressure_win(iy, :, :));
    trace_cut   = squeeze(pressure_win(iy, ix, :));    

    lw = 2;

    figure;
    subplot(3, 1, 1);
    imagesc(time_axis*1e6, 1:Nx, slice_y, [-options.Clip, options.Clip] * max(slice_y(:)));
    xlim(time_axis([1, end]) * 1e6);
    xlabel('Time [us]');
    ylabel('x-position [samples]');
    hold on;
    xline(t_cut*1e6, 'b', 'linewidth', lw);
    yline(ix, 'k', 'linewidth', lw);
    colormap(getBatlow);

    subplot(3, 1, 2);
    imagesc(time_axis_cut*1e6, 1:Nx, slice_y_cut, [-options.Clip, options.Clip] * max(slice_y_cut(:)));
    xlim(time_axis([1, end]) * 1e6);
    xlabel('Time [us]');
    ylabel('x-position [samples]');
    hold on;
    xline(t_cut*1e6, 'b', 'linewidth', lw);
    yline(ix, 'k', 'linewidth', lw);
    colormap(getBatlow);
    
    subplot(3, 1, 3);
    hold on
    plot(time_axis*1e6, trace*1e-3, 'k');
    plot(time_axis_cut*1e6, trace_cut*1e-3, 'r');
    xlim(time_axis([1, end]) * 1e6);
    xlabel('Time [us]');
    ylabel('Pressure [kPa]');
    hold on;
    xline(t_cut*1e6, 'b--', 'linewidth', lw);
    legend({'Original Trace', 'Trimmed Trace'})

    drawnow
end

disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');

end
