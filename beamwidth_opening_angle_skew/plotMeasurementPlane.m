function plotMeasurementPlane(pressure, x_pos, y_pos)

Ny = size(pressure, 1);
Nx = size(pressure, 2);

ssp   = sum(pressure.^2, 3);
max_p = max(pressure, [], 3);

figure; 
subplot(4, 1, 1);
imagesc(x_pos, y_pos, ssp);
axis image;
colorbar;
colormap(getBatlow);
title('Pressure Squared Integral');

subplot(4, 1, 2);
contourf(x_pos, y_pos, ssp/max(ssp(:)), 0.1:0.1:1);
axis image;
colorbar

subplot(4, 1, 3);
imagesc(x_pos, y_pos, max_p);
axis image;
colorbar;
colormap(getBatlow);
title('Maximum Pressure');

subplot(4, 1, 4);
contourf(x_pos, y_pos, max_p/max(max_p(:)), 0.1:0.1:1);
axis image;
colorbar

drawnow

end