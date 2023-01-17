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
