function pressure_win = spatialWindowMeasurementPlane(pressure, border)
%SPATIALWINDOWMEASUREMENTPLANE Window measurement plane in space for every time step
%
% DESCRIPTION:
%     spatialWindowMeasurementPlane applies a 2D spatial window at every
%     time step to an array containing the measured pressure over a plane.
%     The window is defined by its border width, which is the number of
%     samples not equal to one at the edge of the 2D windowing domian. The
%     shape of the window in the border region has is formed using a Tukey
%     window.
%
% USAGE:
%     pressure_win = spatialWindowMeasurementPlane(pressure, border)
%
% INPUTS:
%     pressure     - [numeric] pressure array with size (Ny, Nx, Nt) [Pa]
%     border       - [numeric] integer number of samples defining the
%                    width of the border region
%
% OUTPUTS:
%     pressure_win - [numeric] pressure array with size (Ny, Nx, Nt) [Pa]
%
% ABOUT:
%     author       - Morgan Roberts
%     date         - 2/12/2022

tic;
disp('Apply spatial window to measurement plane ...');

% Work out sizes of each dimension
Ny = size(pressure, 1);
Nx = size(pressure, 2);
Nt = size(pressure, 3);

% Create Tukey window and extract the ring-up and ring-down regions
win = getWin(border*2+1, 'Tukey', 'Param', 1);
ru  = win(1:border);
rd  = win( (border + 2):end );

% Create new windows, controlling the width of the middle region
win_row = [ru; ones(Ny - (2 * border), 1); rd];
win_col = [ru; ones(Nx - (2 * border), 1); rd];

% Use the outer product to create the 2D window
win_2d = win_row * win_col';

% Replicate the window for all time steps
win_3d = repmat(win_2d, [1, 1, Nt]);

% Apply the window
pressure_win = pressure .* win_3d;

disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');

end