% FILENAME:
%     lateral_opening_angle_analysis.m
%
% DESCRIPTION:
%     Script to compute the lateral opening angle of N=48 UST elements.
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

tic;
disp('Loading data ...');
filename       = 'isolated_reprojected_beam_data_probes_AEF';
input_filename = [scan_data_folder, filesep, filename, '.mat'];
load(input_filename, 'beam_data', 'freqs', 'c_water', 'dx', 'Ny', 'Nx', 'Fs', 'dt', 'z_proj', 'z_targ', 'Nfreqs');
disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');

% Calculate x and y position vectors
x_pos = 0:dx:(Nx-1) * dx;
x_pos = x_pos - x_pos( round(Nx / 2) );
y_pos = 0:dx:(Ny-1) * dx;
y_pos = y_pos - y_pos( round(Ny / 2) );
%%

% Setup loop parameters
Nel = size(beam_data, 5);
Nz  = size(beam_data, 4);

% Threshold for defining the edge of the beam
thresh = -6;

% Setup arc interpolation
x_data = repmat( x_pos(:), [1, Nz] );
z_data = repmat( z_targ, [Nx, 1] );
r_every = 1;
Narc   = floor( (z_targ(end) - z_targ(1)) / (r_every*dx) );  
r_q    = z_targ(1) + (0:(r_every*dx):(Narc - 1) * (r_every*dx));
theta  = -60:60;
Ntheta = length(theta);
arc_x  = r_q' * sind(theta);
arc_z  = r_q' * cosd(theta);  % indexed (Nr, Ntheta)   

% Setup storage
dir_resp      = zeros(Ntheta, Nfreqs, Nel);
opening_angle = zeros(1, Nel);

ExtraPlot = 1;

% Loop over each element
for edx = 1:Nel

    % Storage for (Nx, Nf) lateral response spectra, at Nz separate z-planes
    lat_resp = zeros(Nx, Nfreqs, Nz);

    % Extract a lateral response slice for each z-plane
    for zdx = 1:Nz
        % Extract the data for this z-plane: has dims (Ny, Nx, Nfreqs)
        volume_data = squeeze( beam_data(:, :, :, zdx, edx) );
    
        % Find the centroid of the -6dB bounding box in 3D (not frequency
        % restricted)
        [~, Y, ~] = locateThresholdRegion(volume_data, thresh, Round=true);

        % Extract and store the lateral response: has dims (Nx, Nfreqs)
        lat_resp(:,:,zdx) = volume_data(Y.ic, :, :);
    end

    % For each frequency convert the lateral response (cartesian) to a
    % directional response (polar)
    for fdx = 1:Nfreqs
        % Extract lateral resonses for this frequency across all z-planes
        p_data = squeeze( lat_resp(:,fdx,:) );
    
        % Interpolate pressure field along arcs [indexed as (Narc, Ntheta)]
        arc_p = interp2( z_data, x_data, p_data, arc_z, arc_x, 'cubic' );
    
        % Normalise the amplitude along each arc
        arc_p_norm = arc_p ./ repmat( max(arc_p, [], 2, 'omitnan'), [1, Ntheta] );
        
        % Set the amplitude to equal the max amplitude at z_max (z_index=5)
        max_val    = max( squeeze(lat_resp(:,fdx,Nz)), [], 'all' );
        arc_p_norm = arc_p_norm * max_val;
    
        % Calculate the mean response over all arcs for each angle
        arc_p_mean = mean(arc_p_norm, 1, 'omitnan');
    
        % Save directional response for the current frequency and element
        dir_resp(:, fdx, edx) = arc_p_mean;
    
        if ExtraPlot && fdx == 101 && edx == 1
    
            % Interpolated Pressure field for plotting
            p_q = interp2( z_data, x_data, p_data, r_q, x_pos', 'cubic' );
    
            figure;
            subplot(1, 2, 1);
            hold on
            h1 = plot3( 1e3*z_data, 1e3*x_data, p_data, 'r', 'Marker', 'x' );
            h2 = plot3( 1e3*arc_z', 1e3*arc_x', arc_p', 'k', 'Marker', '.' );
            h3 = surf( 1e3*r_q, 1e3*x_pos, p_q, 'EdgeColor', 'none');
            view([27,  15]);
            legend([h1(1), h2(1), h3], {'Measurement points', 'Interpolated Arcs', 'Pressure Field'}, 'location', 'northoutside')
            xlabel('z-position [mm]');
            ylabel('x-position [mm]');
            zlabel('Pressure [Pa]');
        
            subplot(1, 2, 2);
            hold on
            plot(theta, arc_p_norm', 'k.');
            plot(theta, arc_p_mean, 'b-')
            xlabel('Theta [deg]');
            ylabel('Normalised Directional Response');
            set(gcf, 'Position', [316 510 1232 420]);
            legend({'Normalised Interpolated Arcs', 'Mean Response'}, 'location', 'northoutside')

            drawnow

        end
    end

    % Calculate the lateral opening angle
    [F, ~, ~]          = locateThresholdRegion( permute( dir_resp(:,:,edx), [1, 2, 3] ), thresh );
    dir_resp_fc        = dir_resp(:,F.ic,edx);
    opening_angle(edx) = fwhm( dir_resp_fc, theta(2) - theta(1), 0 );

end

% Remove NaN values introduced by interp2 attempting to extrapolate
theta_mask = ~isnan( dir_resp(:,1,1) );
theta      = theta(theta_mask);
Ntheta     = length(theta);
dir_resp   = dir_resp(theta_mask, :, :);

% Calculate the mean lateral opening angle and the standard deviation
mean_oa = mean(opening_angle);
std_oa  = std(opening_angle);
disp(['Lateral opening angle mean = ', num2str(mean_oa, 3),  ' deg, std = ', num2str(std_oa, 3), ' deg']);

% Calculate the mean directional response, convert to dB and compute
% standard deviation
mean_resp = mean(dir_resp, 3, 'omitnan');
max_val   = max(mean_resp, [], 'all', 'omitnan');
mean_dB   = 20 * log10(mean_resp / max_val);
std_resp  = std(dir_resp, 0, 3);
std_pc    = 1e2 * std_resp / max_val;

%%
% ------------------------------------------------------------------------
% Plot the mean lateral opening angle and the standard deviation
plot_thresh = -20;

fig = figure;
subplot(1, 2, 1);
hold on;
imagesc(theta, 1e-6*freqs, mean_dB', [plot_thresh, 0]);
c = colorbar;
set(gca, 'YDir', 'normal');
colormap(getBatlow);
axis square
ylabel('Frequency [MHz]');
xlabel('Theta [deg]');
xlim( theta([1, end]) );
ylim( 1e-6 * freqs([1, end]) );
ylabel(c, 'Pressure [dB]');

subplot(1, 2, 2);
imagesc(theta, 1e-6*freqs, std_pc');
c = colorbar;
set(gca, 'YDir', 'normal');
colormap(getBatlow);
axis square
ylabel('Frequency [MHz]');
xlabel('Theta [deg]');
ylabel(c, 'Standard Deviation [%]');
xlim( theta([1, end]) );
ylim( 1e-6 * freqs([1, end]) );

set(fig, 'Position', [404 278 932 602]);

% ------------------------------------------------------------------------
% Plot the beamwidth histogram

face_colour = [175, 238, 238]/255;

figure;
bin_edges = 48:2:64;
h2 = histogram(opening_angle, bin_edges);
set(h2,'facecolor',face_colour);
xlabel('Opening Angle [deg]');
ylabel('Counts');
axis square
xlim(bin_edges([1, end]))



