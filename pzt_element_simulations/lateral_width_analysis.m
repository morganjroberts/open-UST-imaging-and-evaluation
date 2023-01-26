
close all
clearvars

[~, data_dir]  = getRepoDataPath();
data_folder    = [data_dir, filesep, 'field_scans'];
filename       = 'pzt_element_AFP_amplitude_fields_lateral';
input_filename = [data_folder, filesep, filename, '.mat'];
load(input_filename);


%%


Nx    = kgrid.Nx;
Ny    = kgrid.Ny;
Nz    = kgrid.Nz - z_offset + 1;
z_pos = 0:dx:(Nz - 1) * dx;
x_pos = kgrid.x_vec;

Nfreqs = length(freqs);
Nl     = length(l_width);

% Create query points along an arc centred on the source
r     = z_pos(end);
oa_crit = -90:90;
delta_theta = oa_crit(2) - oa_crit(1);
x_arc = r * sind(oa_crit);
z_arc = r * cosd(oa_crit);

opening_angle = zeros(Nl, Nfreqs);

for ldx = 1:Nl
    for fdx = 1:Nfreqs    
        field = xz_amp_field( :, :, fdx, ldx);

        % Interpolate pressure field along arcs [indexed as (Narc, Ntheta)]
        arc_amp = interp2( z_pos, x_pos, field, z_arc, x_arc, 'cubic' );

        try
            opening_angle(ldx, fdx) = fwhm(arc_amp, delta_theta);
        catch
            opening_angle(ldx, fdx) = 180;
        end
% 
%         figure;
%         subplot(1, 2, 1);
%         surf(z_pos, x_pos, field, 'edgecolor', 'none');
%         hold on
%         plot3(z_arc, x_arc, arc_amp)
% 
%         subplot(1, 2, 2);
%         plot(theta,arc_amp)
%         keyboard



    end
end

% Plot an example of the interpolation process
figure;
subplot(1, 2, 1);
surf(z_pos, x_pos, field, 'edgecolor', 'none');
hold on
plot3(z_arc, x_arc, arc_amp)
subplot(1, 2, 2);
plot(oa_crit,arc_amp)

% plot the final result
figure;
% h = imagesc(freqs*1e-6, l_width*1e3, opening_angle);
contourf(l_width*1e3, freqs*1e-6, opening_angle', [30:10:180]);
c = colorbar;
clim([40, 180])
c.Location = 'northoutside';
colormap(getBatlow);
xlabel('Lateral Width [mm]');
ylabel('Frequency [MHz]');
ylabel(c, 'Opening Angle [deg]');
ylim([0.8, 2]);
axis square


[~,fdx] = findClosest(freqs, 1.2e6);
[~,ldx] = findClosest(l_width, 1e-3);
oa_comp_val = opening_angle(ldx, fdx);
disp(['Predicted opening_angle f = ', ...
    num2str(freqs(fdx)*1e-6), 'MHz: ', num2str(oa_comp_val), '  deg']);

% Calculate mean over frequencies
mean_oa = mean(opening_angle, 2);
figure;
plot(l_width*1e3, mean_oa, 'k-', 'linewidth', 1.5);
% axis square
xlim( 1e3*l_width([1, end]) );
xlabel('Lateral Width [mm]');
ylabel('Mean Opening Angle [deg]')
set(gcf, 'Position', [336 1027 389 194]);

r_roi = 0:110;
r_array = 110;
oa_crit = 2 * asind(r_roi ./ r_array);

D_roi_max = interp1(oa_crit, r_roi, opening_angle(ldx,:)) * 2;


figure;
plot(freqs, D_roi_max)
D_roi_crit = min(D_roi_max);


% disp(['Max object diameter for opening angle of ', ...
%      num2str(oa_comp_val), ' deg: ', num2str(D_roi_max), ' mm']);

disp(['Max object diameter for lateral width of ', ...
     num2str(l_width(ldx)*1e3), ' mm: ', num2str(D_roi_crit), ' mm']);


D_roi_mean = interp1(oa_crit, r_roi, oa_comp_val) * 2;
disp(['Object diameter for opening angle of ', ...
     num2str(oa_comp_val), '  deg: ', num2str(D_roi_mean), ' mm']);