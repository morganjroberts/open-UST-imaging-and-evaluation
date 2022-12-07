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

beam_data = squeeze(beam_data(:,:,:,:,1));

tic;
disp('Saving data ...');
filename        = 'isolated_reprojected_beam_data_ch1';
output_filename = [scan_data_folder, filesep, filename, '.mat'];
keyboard
save(output_filename, 'beam_data', 'freqs', 'c_water', 'dx', 'Ny', 'Nx', 'Fs', 'dt', 'z_proj', 'z_targ', 'Nfreqs');
disp(['Completed in ', num2str(toc), ' s']);
fprintf('\n');