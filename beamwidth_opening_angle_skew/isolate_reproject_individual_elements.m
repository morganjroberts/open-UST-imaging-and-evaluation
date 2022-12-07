% FILENAME:
%     isolate_reproject_individual_elements.m
%
% DESCRIPTION:
%    
%
% APPROACH:
%     
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
filenames       = {'field_scan_probe_A_all_channels_postprocessed', ...
                   'field_scan_probe_E_all_channels_postprocessed', ...
                   'field_scan_probe_F_all_channels_postprocessed'};

input_filename = [scan_data_folder, filesep, filenames{1}, '.mat'];
load(input_filename, 'source_p', 'source_z', 'c_water', 'dx', 'dt', 'z_pos');
loaded_probe = 1;

% Target plane locations and projection distances
z_targ = 70e-3:10e-3:110e-3;
z_proj = source_z + z_targ - abs(z_pos);

% Experiment parameters
Nelement = 48; 
Nx_proj  = 401;
Ny_proj  = 79;
Nz_pos   = length(z_proj);

Nper_mod = 16;
Nmodules = Nelement / Nper_mod;

% Work out how many frequencies to allocate space for
f_min     = 0.4e6;
f_max     = 2.0002e6;
f_range   = f_max - f_min;
Nfreqs    = 201; % ideal number of frequencies
Fs        = 1 / dt;
Fnyq      = Fs / 2;
delta_f   = round(f_range / (Nfreqs - 1));
Nfft_half = round(Fnyq / delta_f);
Nfft      = 2 * Nfft_half;

% Create a test vector to check the frequencies
a           = randn(1, Nfft);
[freqs, ~]  = spect(a, Fs, 'FFTLength', Nfft);

% Work out the actual number of frequencies
freqs(freqs < f_min) = 0;
freqs(freqs > f_max) = 0;
i_f_start            = find(freqs, 1, 'first');
i_f_end              = find(freqs, 1, 'last');
freqs                = freqs(i_f_start:i_f_end);
Nfreqs               = length(freqs);
delta_f              = freqs(2) - freqs(1);

% Check the widest opening angle that is contained by the proposed grid
% length
Nz_max     = round(z_proj(end) / dx);
theta_max  = 2 * atand(Nx_proj / (2 * Nz_max));
disp(['Max opening angle at z = ', num2str(z_proj(end)*1e3), ' mm is ', ...
    num2str(theta_max), ' deg']);

% Initialise storage array, check it leaves enough memory headroom
beam_data = zeros(Ny_proj, Nx_proj, Nfreqs, Nz_pos, Nelement, 'single');
mem       = whos('beam_data');
disp(['Storage required: ', num2str(mem.bytes), ' bytes']);
fprintf('\n');

% Setup initial estimates for element centroids (visual assessment)
plotMeasurementPlane(source_p)
ix0_ch1   = 199;
pitch     = 2.54e-3;
pitch_pts = pitch / dx;
ix0s      = round(ix0_ch1 - (0:pitch_pts:((Nper_mod - 1) * pitch_pts)));
iy0       = 31;

% Parameters for extractSourceElement
border_param = [10,80];
phys_size    = [0.8e-3, 9.8e-3];
ups          = 7;
mask_size    = [0.00235, 0.01135];


as_options  = {'GridExpansion', 0, 'Plot', 0, 'Reverse', 0, 'DataCast', 'single', 'DataRecast', 0, 'FFTLength', 4096};

for idx = 1:Nelement
    [channel, probe_id] = ind2sub([Nper_mod, Nmodules], idx);
    disp(['Probe ', num2str(probe_id), ', channel ', num2str(channel)]);

    % Check the correct dataset is loaded. If not, load it and update flag.
    if loaded_probe ~= probe_id
        input_filename = [scan_data_folder, filesep, filenames{probe_id}, '.mat'];
        disp('New probe dataset required, loading in progress ...');
        load(input_filename, 'source_p', 'c_water');
        loaded_probe = probe_id;
    end

    % Extract the source pressure from a single channel
    extract_p = extractSourceElement(source_p, dx, [iy0, ix0s(channel)], ...
        Border=border_param, UpSample=ups, ExtraPlot=false, ...
        PhysicalSize=phys_size, MaskSize=mask_size);

    % Create a new domain and centre the element source pattern within it
    centred_p = centreElementInNewGrid(extract_p, Nx_proj, Ny_proj, ExtraPlot=false);

    % Perform the reprojection
    [~, proj_p] = angularSpectrum(centred_p, dx, dt, z_proj, c_water, as_options{:});
    % (y_ind, x_ind, t_ind, z_ind)

    % Compute and store amplitude spectra
    [~, as_meas] = spect(proj_p, Fs, 'FFTLength', Nfft, 'Dim', 3);
    beam_data(:,:,:,:,idx) = as_meas(:,:,i_f_start:i_f_end,:);

    % Append to file
    Ny = Ny_proj;
    Nx = Nx_proj;
    filename = 'isolated_reprojected_beam_data_probes_AEF';
    output_filename = [scan_data_folder, filesep, filename, '.mat'];
    if ~exist(output_filename, 'file')
        save(output_filename, 'beam_data', 'freqs', 'c_water', 'dx', 'dt', 'Ny', 'Nx', 'Fs', 'Nfreqs', 'z_proj', 'z_targ', '-v7.3');
    else 
        save(output_filename, 'beam_data', 'freqs', 'c_water', 'dx', 'dt', 'Ny', 'Nx', 'Fs', 'Nfreqs', 'z_proj', 'z_targ', '-append');
    end
end



function centred_p = centreElementInNewGrid(source_p, Nx_new, Ny_new, options)

arguments
    source_p
    Nx_new
    Ny_new
    options.ExtraPlot = true;
end

Nt = size(source_p, 3);

% Find the location of the non-zero region of source_p
ssp = sum(source_p.^2, 3);
iL  = find(sum(ssp, 1), 1, 'first');
iR  = find(sum(ssp, 1), 1, 'last');
iT  = find(sum(ssp, 2), 1, 'first');
iB  = find(sum(ssp, 2), 1, 'last');

% size of the non-zero region
Nx_extr = iR - iL + 1;
Ny_extr = iB - iT + 1;

% create new computational grid
centred_p = zeros(Ny_new, Nx_new, Nt);

% location within the new grid to place the extracted region 
x_start = ceil(Nx_new / 2) - round(Nx_extr / 2);
x_end   = x_start + Nx_extr - 1;
y_start = ceil(Ny_new / 2) - round(Ny_extr / 2);
y_end   = y_start + Ny_extr - 1;

% Copy the extracted pressure into the centre of the new grid
centred_p( y_start:y_end,  x_start:x_end, : ) = source_p( iT:iB, iL:iR, : );

if options.ExtraPlot
    ssp_new = sum(centred_p.^2, 3);

    figure;
    hold on
    imagesc(ssp_new);
    axis image
    h = plot(round(Nx_new/2), round(Ny_new/2), 'kx');
    c = colorbar;
    colormap(getBatlow);
    legend([h], {'Centre of Domain'});
    xlabel('x-position [mm]');
    ylabel('y-position [mm]');
    ylabel(c, 'Pressure squared integral');
end


end


% setup linear indexing for 16 x 3 to load proper file
% setup a vector of projected z-positions based on 4th gen
% for each element, 
% - load x-location and y-location
% - work out whether the currently loaded file is valid
% - extractSourceElement
% - centre the element in the middle of the array
% - reproject to the z-positions
% - compute amplitude spectrum (restricted to only essential freqs for
% memory)
% - store as (Ny, Nx, Nf, Nz, Nel)



