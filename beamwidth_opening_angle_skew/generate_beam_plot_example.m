% FILENAME:
%     isolation_reprojection_mask_size_validation.m
%
% DESCRIPTION:
%     It is required to assess the beam shapes of N=48 elements. Instead of
%     perofrming field scans for multiple elements, entire modules were
%     scanned with all elements transmitting, to save time and memory. To
%     estimate the beam shapes of individual elements from a scan with
%     multiple elements being driven, the source pressure distribution of
%     each element must be isolated at the source plane, using a mask.
%     Since the mask design affects the field shape, the mask design must
%     be validated. As a comparison, a field scan was performed where only
%     a single element is transmitting (probe A channel 8). Also in this
%     script, it is shown that grid expansion can be set to zero for the
%     angular spectrum reprojections, to save time and memory, with
%     negligible change in error.
%
% APPROACH:
%     Pre-requisite: Project measurement plane back to source plane
%     (completed in backproject_to_source_plane.m)
%     1. Create an array of different mask sizes. 
%     2. Use the function extractSourceElement.m to extract the source
%     pressure distribution of channel 8 only, for each of the mask sizes.
%     3. For each separate extracted pressure, re-project to the validation
%     plane (z = 80.15mm)
%     4. Calculate the normalised amplitude at 1.2MHz for this re-projected plane
%     5. Also calculate the normalised amplitude at 1.2MHz for the validation scan
%     6. Store the sum of absolute differences between (4) and (5).
%     7. Repeat for all mask sizes. The mask sizes producing the lowest
%     error is optimal, and will be used for the isolation/reprojection
%     pipeline to estimate beamshapes of individual elements.
% 
% INPUT DATA FILENAMES:
%     <data-dir>\field_scans\_______.mat      
%
% ABOUT:
%     author:      - Morgan Roberts
%     last update: - 6/12/22

close all
clearvars

[~, data_dir]    = getRepoDataPath();
scan_data_folder = [data_dir, filesep, 'field_scans'];

% -------------------------------------------------------------------------
% Load the datasets and backpropagate to a series of planes
filename       = 'field_scan_probe_A_all_channels_postprocessed';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
load(input_filename, 'source_p', 'source_z', 'c_water', 'dx', 'dt', 'z_pos');

% For UST probes in UCL scan tank, channel 1 has most +ve x-coordinate
% Visual identification of element centre (ch8) after border cutting
plotMeasurementPlane(source_p);
ix0 = 148; 
iy0 = 31; 

% The real source is located at z=0, but the optimal source plane is
% located at z = 1.4 mm
% source_z is the distance from the measurement plane to the optimal source
% plane
z_s = abs(z_pos) - source_z;

% The final plane position, relative to z = 0
z_targ = ceil(110e-3 / dx) * dx;
% Final plane relative to optimal source plane
z_proj_end = z_targ - z_s;
% All plane positions relative to optimal source plane
z_proj = -z_s:dx:z_proj_end;

% Setup angular spectrum parameters
as_options   = {'GridExpansion', 0, 'Plot', 0, 'Reverse', 0, 'DataCast', 'single', 'DataRecast', 1};

% Create range of mask sizes
phys_size = [0.8e-3, 9.8e-3];
mask_size = [2.35e-3, 11.35e-3];

% Extract the source pressure from a single channel
extract_p = extractSourceElement(source_p, dx, [iy0, ix0], ...
    Border=[10,80], UpSample=7, ExtraPlot=false, ...
    PhysicalSize=phys_size, MaskSize=mask_size);

% Create a new domain and centre the element source pattern within it
Nx_proj  = 401;
Ny_proj  = 79;
centred_p = centreElementInNewGrid(extract_p, Nx_proj, Ny_proj, ExtraPlot=false);
    
% Project the source to the measurement plane
max_p = angularSpectrum(centred_p, dx, dt, z_proj, c_water, as_options{:});
    
filename_out    = 'probe_A_channel_8_maximum_pressure_field';
output_filename = [scan_data_folder, filesep, filename_out, '.mat'];
save(output_filename, 'max_p', 'dx', 'source_z', 'z_pos', 'z_proj');

%%



Nx = size(max_p, 2);
Ny = size(max_p, 1);
Nz = size(max_p, 3);

x_pos = 0:dx:(Nx - 1)*dx;
x_pos = x_pos - x_pos( round(Nx/2) );
y_pos = 0:dx:(Ny - 1)*dx;
y_pos = y_pos - y_pos( round(Ny/2) );
z_pos = 0:dx:(Nz - 1)*dx;

yz_plane = squeeze( max_p(:, round(Nx/2), :) );
xz_plane = squeeze( max_p(round(Ny/2), :, :) );
xy_plane = squeeze( max_p(:,:,end) );

a = 100;
% scale data between 0 and 1
yz = (yz_plane - min(yz_plane(:))) ./ max(yz_plane(:));
xz = (xz_plane - min(xz_plane(:))) ./ max(xz_plane(:));
xy = (xy_plane - min(xy_plane(:))) ./ max(xy_plane(:));
% compress
yz_lc = log10(1 + a * yz) ./ log10(1 + a);
xz_lc = log10(1 + a * xz) ./ log10(1 + a);
xy_lc = log10(1 + a * xy) ./ log10(1 + a);

r = 103e-3;
theta  = -60:2:60;
arc_x  = r * sind(theta);
arc_z  = r * cosd(theta);
valid = (arc_z > 70e-3) & (arc_x > min(x_pos)) & (arc_x < max(x_pos));
arc_x = arc_x(valid);
arc_z = arc_z(valid);

i_cut = 4;

figure;
subplot(1, 2, 2);
imagesc(y_pos*1e3, z_pos(i_cut:end)*1e3, yz_lc(:,i_cut:end)');
yline([70:10:110], 'k-.', 'linewidth', 1.5);
c = colorbar(gca, 'eastoutside');
colormap(getColorMap);
axis image
xlabel('y-position [mm]');

subplot(1, 2, 1);
imagesc(x_pos*1e3, z_pos(i_cut:end)*1e3, xz_lc(:,i_cut:end)');
hold on
plot(arc_x*1e3, arc_z*1e3, 'k-', 'linewidth', 1.5);
yline([70:10:110], 'k-.', 'linewidth', 1.5);
colormap(getColorMap);
axis image
ylabel('z-position [mm]');
xlabel('x-position [mm]');

set(gcf, 'Position', [-1862 1197 1735 515]);
%%
% TUFFC2023
% Generate elevational profile and directional response profile.
% These will be added manually to the xz and yz slices in inkscape

figure;
beam_profile = squeeze( yz_plane(:,end) );
beam_profile = beam_profile / max(beam_profile);

fwhm(beam_profile, 1, 1)
xlim([1, Ny-1])
ylim([0.1, 1.05])

arc_p = interp2( z_pos, x_pos, xz_plane, arc_z, arc_x, 'cubic' );
arc_p = arc_p / max(arc_p);
figure;
fwhm(arc_p, 1, 1)
xlim([1, length(arc_z)-1])
ylim([0.1, 1.05])

