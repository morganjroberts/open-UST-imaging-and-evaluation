

% SIMULATIONS TO MEASURE EFFECT OF ELEVATION HEIGHT AND FREQUENCY ON
% BEAMWIDTH
% Already established that lateral width affects the field amplitude, but
% the beamwidht is unaffected (from 0.7mm - 1.4mm)

close all
clearvars

[~, data_dir] = getRepoDataPath();
data_folder   = [data_dir, filesep, 'field_scans'];
        
% Create arrays for the various combinations of element dimensions
e_height = 1e-2;
l_width  = 1e-4 * (7:15);
Nl       = length(l_width);

% create homogeneous lossless medium
temperature = 22;
c           = waterSoundSpeed(temperature);

% calculate grid size
f_max = 2.2e6;
ppw   = 3;
dx    = c / (ppw * f_max);

% create kgrid
Nx       = 495;
Ny       = 90;
Nz       = 250;
z_offset = 15;
kgrid    = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

disp( ['Lx: ', num2str( 1e3* (Nx - 1) * dx ), ' mm'] )
disp( ['Ly: ', num2str( 1e3* (Ny - 1) * dx ), ' mm'] )
disp( ['Lz: ', num2str( 1e3* (Nz - 1 - z_offset) * dx ), ' mm'] )

%     x_size = 110e-3;
%     y_size = 35e-3;
%     z_size = 110e-3;

%%

% off grid source setup
cx       = kgrid.x_vec( round( Nx / 2 ) );
cy       = kgrid.y_vec( round( Ny / 2 ) );
cz       = kgrid.z_vec( z_offset );
position = [cx, cy, cz];
theta    = [0, 0, 0];

% create frequency vector
freqs = 0.8e6:60e3:2.1e6;
Nf    = length(freqs);

% initialise storage
xz_amp_field = zeros(Nx, Nz-z_offset+1, Nf, Nl, 'single');

for fdx = 1:Nf
    for ldx = 1:Nl
    
            % Create new kWaveArray
            karray = kWaveArray('UpsamplingRate', 500);
    
            % Create offgrid source with new element dimensions
            Lx = l_width(ldx);
            Ly = e_height;
            karray.addRectElement(position, Lx, Ly, theta)
            amp_in = karray.getArrayGridWeights(kgrid);
          
            phase_in = 0;
            [amp_out, ~] = acousticFieldPropagatorC(amp_in, phase_in, dx, ...
                           freqs(fdx), c);
    
            xz_amp_field( :, :, fdx, ldx) = amp_out( :, round(Ny/2), z_offset:end);
 
            data = squeeze( amp_out( :, round(Ny/2), z_offset:end) );
            data_db = 20 * log10( data / max( data(:,end) ) );

%             z_pos = kgrid.z_vec(z_offset:end) - kgrid.z_vec(z_offset);
%             x_pos = kgrid.x_vec;
%             figure;
%             imagesc(z_pos, x_pos, data );
%             hold on;
%             colormap(getBatlow)
%             axis image
%             drawnow

%             keyboard


% 
%             figure;
%             plot(z_pos, calculateBeamwidthProfile(data, dx))
%             drawnow

            
    end
end

filename = 'pzt_element_AFP_amplitude_fields_lateral';
output_filename = [data_folder, filesep, filename, '.mat'];
save(output_filename, 'xz_amp_field', 'freqs', 'e_height',...
    'l_width', 'kgrid', 'z_offset', 'dx', '-v7.3');



