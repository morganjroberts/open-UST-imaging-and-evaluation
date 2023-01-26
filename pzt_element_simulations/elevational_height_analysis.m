
close all
clearvars

[~, data_dir]  = getRepoDataPath();
data_folder    = [data_dir, filesep, 'field_scans'];
filename       = 'pzt_element_AFP_amplitude_fields_elevational';
input_filename = [data_folder, filesep, filename, '.mat'];
load(input_filename);


%%

Nz = kgrid.Nz;
Nx = kgrid.Nx;
Ny = kgrid.Ny;

Nz = Nz - z_offset + 1;

z_pos = 0:dx:(Nz - 1) * dx;
x_pos = 0:dx:(Nx - 1) * dx;
x_pos = x_pos - x_pos( round(Nx/2) );
y_pos = 0:dx:(Ny - 1) * dx;
y_pos = y_pos - y_pos( round(Ny/2) );

Nfreqs = length(freqs);
Ne = length(e_height);

% find far field z-position
[~, z_end]   = findClosest(z_pos, 0.11);
[~, z_start] = findClosest(z_pos, 0.0);
i_z = z_start:z_end;

mean_beamwidth = zeros(Ne, Nfreqs);
beamwidth_c = zeros(Ne, Nfreqs);

for edx = 1:Ne
    for fdx = 1:Nfreqs    
        bws = zeros((z_end-z_start+1), 1);
        for zdx = 1:length(i_z)
            ele_profile = squeeze( yz_amp_field(:, i_z(zdx), fdx, edx) );
            
            try
                bws(zdx) = fwhm(ele_profile, dx);
            catch
                bws(zdx) = NaN;
            end
        end

        mean_beamwidth(edx, fdx) = mean(bws);
        beamwidth_c(edx, fdx) = bws(end);
% 
%         figure;
%         plot(bws);
%         keyboard


    end
end

% [~,fdx] = findClosest(freqs, 1.2e6);
% data = beamwidth(:,:,fdx);
% figure;
% h = imagesc(l_width*1e3, e_height*1e3, data*1e3);
% c = colorbar;
% colormap(getBatlow);
% set(h, 'AlphaData', ~isnan(data));
% ylabel('Elevation Height [mm]');
% xlabel('Lateral Width [mm]');
% ylabel(c, 'Beamwidth [mm]');

% Plot effect of frequency and elevation height on beamwidth

% work out optimal e_height for every f
[~,I] = min(mean_beamwidth, [], 1);
e_opt = e_height(I);

figure;
plot(e_height*1e3, mean(mean_beamwidth, 2)*1e3, 'k-', 'linewidth', 1.5 );
xlim( 1e3*e_height([1, end]) );
xlabel('Elevation Height [mm]');
ylabel('Mean Beamwidth [mm]')
set(gcf, 'Position', [336 1034 378 187]);

figure;
hold on
% imagesc(e_height*1e3, freqs*1e-6, mean_beamwidth'*1e3);
contourf( e_height*1e3, freqs*1e-6, mean_beamwidth'*1e3, 7:0.5:18);
clim([8,18]);
h = plot(e_opt*1e3, freqs*1e-6, 'rx');
h.MarkerEdgeColor = 'w';
c = colorbar;
c.Location = 'northoutside';
colormap(getBatlow);
xlabel('Elevation Height [mm]');
ylabel('Frequency [MHz]');
ylabel(c, 'Beamwidth [mm]');
ylim([0.8, 2]);
axis square
box on

figure;
imagesc(freqs*1e-6, e_height*1e3, beamwidth_c*1e3);
% h = contourf(freqs*1e-6, e_height*1e3, beamwidth_c*1e3);
c = colorbar;
colormap(getBatlow);
ylabel('Elevation Height [mm]');
xlabel('Frequency [MHz]');
ylabel(c, 'Beamwidth [mm]');
axis square

[~,fdx] = findClosest(freqs, 1.2e6);
[~,edx] = findClosest(e_height, 10e-3);
bw_comp_val = beamwidth_c(edx, fdx);
disp(['Predicted beamwidth at z = 110 mm, f = ', ...
    num2str(freqs(fdx)*1e-6), 'MHz: ', num2str(bw_comp_val*1e3), '  mm']);

figure;
plot(freqs*1e-6, mean_beamwidth(edx+1,:)*1e3)
