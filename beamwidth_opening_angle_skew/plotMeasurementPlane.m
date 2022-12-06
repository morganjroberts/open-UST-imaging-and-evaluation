function plotMeasurementPlane(pressure, options)

arguments
    pressure
    options.Xpos = 1:size(pressure, 2);
    options.Ypos = 1:size(pressure, 1);
end

Ny = size(pressure, 1);
Nx = size(pressure, 2);

ssp   = sum(pressure.^2, 3);
max_p = max(pressure, [], 3);

figure; 
subplot(4, 1, 1);
imagesc(options.Xpos, options.Ypos, ssp);
axis image;
colorbar;
colormap(getBatlow);
title('Pressure Squared Integral');
set(gca, 'YDir', 'normal');

subplot(4, 1, 2);
contourf(options.Xpos, options.Ypos, ssp/max(ssp(:)), 0.1:0.1:1);
axis image;
colorbar
set(gca, 'YDir', 'normal');

subplot(4, 1, 3);
imagesc(options.Xpos, options.Ypos, max_p);
axis image;
colorbar;
colormap(getBatlow);
title('Maximum Pressure');
set(gca, 'YDir', 'normal');

subplot(4, 1, 4);
contourf(options.Xpos, options.Ypos, max_p/max(max_p(:)), 0.1:0.1:1);
axis image;
colorbar
set(gca, 'YDir', 'normal');

drawnow

end