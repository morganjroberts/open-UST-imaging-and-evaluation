function theta = calculateBulkLateralSkew(data_dir, filename, options)

arguments
    data_dir
    filename
    options.ExtraPlot = true;
end

input_filename = [data_dir, filesep, filename, '.mat'];
load(input_filename, 'source_p', 'time_axis', 'c_water', 'x_pos', 'Ny');

% estimate where the elements are
ssp      = sum( source_p.^2, 3 );
ssp_line = ssp( round(Ny/2), : );
ssp_line = ssp_line / max(ssp_line);
sspl     = ssp_line;
thresh   = 0.4;
sspl(sspl < thresh) = 0;

% find the indexes corresponding to each element
[pks,locs] = findpeaks(sspl, 'NPeaks', 16);

% extract one pressure trace from each element
extr_p = squeeze( source_p(round(Ny/2), locs, :) );

% find the location of the maximum pressure in each trace
[~, I] = max(extr_p, [], 2);

% calculate the relative pressure-peak arrival time of each element
t_arr = time_axis(I);
t_rel = t_arr - t_arr(1);

% convert relative arrival time to relative z-position
z_rel = c_water * t_rel;

% calculate the x-cordinate of each element
x_el = x_pos(locs);

% fit a line
P = polyfit(x_el, z_rel, 1);
z_fit = polyval(P, x_el);

% calculate bulk skew angle from the gradient
theta = atand( P(1) );

% this method uses the convention that a positive bulk skew produces a
% negative bulk skew in the beam axis (ve correlation between x and z)
% To remove the bulk skew later, subtract it.
% correct for this:
theta = -theta;

% plot the process
if options.ExtraPlot
    figure;
    subplot(2, 2, 3);
    hold on;
    plot(1e3*x_pos, ssp_line, 'k');
    plot(1e3*x_pos(locs), pks*max(ssp_line), 'rx');
    xlabel('x-position [mm]');
    
    subplot(2, 2, [1, 2]);
    imagesc(1e6*time_axis, 1:16, 1e-3*extr_p);
    hold on;
    h = plot(1e6*time_axis(I), 1:16, 'r');
    h.Marker = 'x';
    c = colorbar;
    colormap(getBatlow)
    xlim( 1e6*time_axis([2000, 2150]) );
    xlabel('Time [us]');
    ylabel('Element index');
    ylabel(c, 'Pressure [kPa]');
    
    subplot(2, 2, 4);
    hold on
    plot(x_el*1e3, z_rel*1e3, 'rx');
    plot(x_el*1e3, z_fit*1e3, 'b');
    xlabel('x-position [mm]');
    ylabel('z-position [mm]');
end

end