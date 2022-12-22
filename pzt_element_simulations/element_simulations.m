% set up kwave simulation offgrid source
% multiple source sizes

% send in all frequencies
% store time series at every point
% calculate amplitude spectrum at every point

% use input signal FFT to scale amplitude field field

close all
clearvars

[~, data_dir] = getRepoDataPath();
data_folder   = [data_dir, filesep, 'field_scans'];
        
% Create arrays for the various combinations of element dimensions
e_height = 1e-3 * (7:14);
l_width  = 1e-4 * (7:14);
Ne       = length(e_height);
Nl       = length(l_width);

% create homogeneous lossless medium
temperature        = 22;
medium.sound_speed = waterSoundSpeed(temperature);

% calculate grid size
f_max = 2.2e6;
ppw   = 3;
dx    = medium.sound_speed / (ppw * f_max);

% create kgrid
Nx    = 648;
Ny    = 120;
Nz    = 540;
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% Create sensor mask to save pressure on planes through element centroid
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(round(Nx/2),:,:) = 1;
sensor.mask(:,round(Ny/2),:) = 1;

% Create vectorised mask to extract xz and yz planes separately
xz_mask = zeros(Nx, Ny, Nz);
yz_mask = zeros(Nx, Ny, Nz);
xz_mask( :, round(Ny/2), : ) = 1;
yz_mask( round(Nx/2), :, : ) = 1;
xz_mask = xz_mask(:);
yz_mask = yz_mask(:);
% Find linear indices of the sensor points in each plane
xzI = find(xz_mask);
yzI = find(yz_mask);
I   = find(sensor.mask);
% Find where the plane data will be located within the sensor data
xz_loc  = ismember(I, xzI);
yz_loc  = ismember(I, yzI);

% Create time array
CFL           = 0.3;
t_end         = 1.8 * 110e-3 /  medium.sound_speed;
kgrid.t_array = makeTime(kgrid, medium.sound_speed, CFL, t_end);

% Create excitation signal, filtered so it is supported on grid
source_f   = 1.2e6;
Fs         = 1 / kgrid.dt;
source_sig = toneBurst(Fs, source_f, 2, 'Plot', true, 'SignalLength', kgrid.Nt);
source_sig = filterTimeSeries(kgrid, medium, source_sig, ...
    'PlotSignals', true, 'PlotSpectrums', true, 'StopBandAtten', 80, ...
    'TransitionWidth', 0.01);

% off grid source setup
z_offset = 30;
cx       = kgrid.x_vec( round( Nx / 2 ) );
cy       = kgrid.y_vec( round( Ny / 2 ) );
cz       = kgrid.z_vec( z_offset );
position = [cx, cy, cz];
theta    = [0, 0, 0];

% Set PML size
pml_size = 10;

% ------------------------------------------------------------------------
% Loop over element sizes, compute amplitude field, and store

% Work out how many frequencies to allocate space for
f_min     = 0.4e6;
f_max     = 2.1e6;
f_range   = f_max - f_min;
Nfreqs    = 201; % ideal number of frequencies
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
freqs_ext            = freqs(i_f_start:i_f_end);
Nfreqs               = length(freqs_ext);
delta_f              = freqs_ext(2) - freqs_ext(1);

% Indexing vector to select only some frequencies to save (for memory)
Nf_save  = 50;
i_f_step = ceil(Nfreqs / Nf_save);
i_f_vec  = i_f_start:i_f_step:i_f_end;
freqs    = freqs(i_f_vec);

% Adjust grid dimensions with removed PML and z_offset
i_x = (pml_size + 1):(Nx - pml_size); 
i_y = (pml_size + 1):(Ny - pml_size); 
i_z = z_offset:(Nz - pml_size);
Nx2 = length(i_x);
Ny2 = length(i_y);
Nz2 = length(i_z);

% initialise storage
xz_amp_field = zeros(Nx2, Nz2, Nf_save, Ne, Nl, 'single');
yz_amp_field = zeros(Ny2, Nz2, Nf_save, Ne, Nl, 'single');

for edx = 1:Ne
    for ldx = 1:Nl

        % Create new kWaveArray
        karray = kWaveArray('UpsamplingRate', 500);

        % Create offgrid source with new element dimensions
        Lx = l_width(ldx);
        Ly = e_height(edx);
        karray.addRectElement(position, Lx, Ly, theta)
        source.p_mask = karray.getArrayBinaryMask(kgrid);

        % Assign source pressure to the off grid sources (BLI)
        source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
        
        % Run simulation
        options     = {'PMLSize', pml_size};
        sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, options{:});
        
        % Distribute the 3D sensor data to the separate 2D planes
        xz_data = sensor_data( xz_loc , :);
        yz_data = sensor_data( yz_loc , :);
        
        % get some memory back
        clear sensor_data

        % Reshape the sensor data to 2D planes
        xz_data = reshape(xz_data, [Nx, Nz, kgrid.Nt]);
        yz_data = reshape(yz_data, [Ny, Nz, kgrid.Nt]);
                
        % Compute the amplitude fields
        [~, xz_as] = spect(xz_data, Fs, 'Dim', 3, 'FFTLength', Nfft);
        [~, yz_as] = spect(yz_data, Fs, 'Dim', 3, 'FFTLength', Nfft);
        
        % Store the amplitude fields with PML and z_offset removed
        xz_amp_field( :, :, :, edx, ldx) = xz_as ( i_x, i_z, i_f_vec );
        yz_amp_field( :, :, :, edx, ldx) = yz_as ( i_y, i_z, i_f_vec );

        % get some memory back
        clear xz_as yz_as

    end
end

filename = 'pzt_element_simulated_amplitude_fields';
output_filename = [data_folder, filesep, filename, '.mat'];
save(output_filename, 'xz_amp_field', 'yz_amp_field', 'freqs', 'e_height',...
    'l_width', 'kgrid', 'Nx2', 'Ny2', 'Nz2', 'dx', 'medium', 'source', '-v7.3');







