% FILENAME:
%     generate_phantom_true_sound_speed.m
% 
% DESCRIPTON:
%     Use dimensions from the CAD model and estimated true sound speed
%     values to build a numerical version of the UST phantom, for
%     comparison with image reconstruction estimates.
%
% REFERENCE DATA FILENAMES:
%     <repo-dir>\phantom_design\UST_phantom.f3d                          CAD model for phantom
%     <repo-dir>\phantom_design\phantom_true_sound_speed_measurement.m   Analysis to estimate the true sound speeds of the liquids used to build the phantom
%     
% OUTPUT DIRECTORY:
%     <repo-dir>\phantom_design
%
% DEPENDENCIES:
%     k-Wave   https://github.com/ucl-bug/k-wave
%
% ABOUT:
%     author:      - Morgan Roberts
%     last update: - 17/11/22

close all
clearvars

% Set grid dimensions
Lx = 152;               % Physical grid length [mm]
dx = 0.25;              % Step size [mm]
Nx = round(Lx/dx) + 1;  % Grid length [pts]

% Index corresponding to centre of grid
pos0 = round(Nx/2);

% Create grid vector
grid_x = 0:dx:(Nx-1)*dx;
grid_x = grid_x - grid_x(pos0);

% Create binary mask features for the different phantom layers (dimensions from CAD model)
oil   = logical(makeDisc(Nx, Nx, 0,                   0,                    round(49/dx)));
fibro = logical(makeDisc(Nx, Nx, round(pos0+6/dx),    round(pos0+0.365/dx), round((69/2)/dx)));
tumA  = logical(makeDisc(Nx, Nx, round(pos0+30.3/dx), round(pos0+0.365/dx), round(4.5/dx)));
tumB  = logical(makeDisc(Nx, Nx, round(pos0-6/dx),    round(pos0+0.365/dx), round(4.5/dx)));
tumC  = logical(makeDisc(Nx, Nx, round(pos0-18.3/dx), round(pos0+0.365/dx), round(2.5/dx)));
tumD  = logical(makeDisc(Nx, Nx, round(pos0-6/dx),    round(pos0+22.51/dx), round(2.5/dx)));
tumE  = logical(makeDisc(Nx, Nx, round(pos0-6/dx),    round(pos0-21.78/dx), round(4.5/dx)));
tumF  = logical(makeDisc(Nx, Nx, round(pos0+16/dx),   round(pos0+22.51/dx), round(4.5/dx)));
tumG  = logical(makeDisc(Nx, Nx, round(pos0+16/dx),   round(pos0-21.78/dx), round(2.5/dx)));

% Assign sound speed values to the different regions, starting with the
% background, moving to the foreground
c0             = waterSoundSpeed(20);
phantom        = c0 * ones(Nx, Nx);
phantom(oil)   = 1457;
phantom(fibro) = 1518; 
phantom(tumA)  = 1526;   % tumour sound speed group 1
phantom(tumG)  = 1526;   % tumour sound speed group 1
phantom(tumB)  = 1535;   % tumour sound speed group 2
phantom(tumC)  = 1535;   % tumour sound speed group 2
phantom(tumF)  = 1535;   % tumour sound speed group 2
phantom(tumD)  = 1544;   % tumour sound speed group 3  
phantom(tumE)  = 1544;   % tumour sound speed group 3

% Correct the orientation of the phantom to match experimental setup
phantom = rot90(phantom', 2);

% Plot figure
fig = figure;
imagesc(grid_x, grid_x, phantom, [1420, 1580]);
colormap(getColorMap);
c = colorbar;
xlabel('x-position [mm]');
ylabel('y-position [mm]');
ylabel(c, 'Sound Speed [m/s]')
axis image

% Save figure
set(fig,'renderer','Painters');
figure_filename = [pwd, '\phantom_true_sound_speed'];
print(fig, figure_filename, '-depsc2');
print(fig, figure_filename, '-dsvg');

% Save the phantom
save('phantom_true_sound_speed.mat', 'phantom', 'dx');
