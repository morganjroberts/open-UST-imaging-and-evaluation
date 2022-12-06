function mask = createElementMask(Nx, Ny, dx, sourceW, sourceL, ix0, iy0)


% check that dx is a factor of sourceW and sourceL
if rem(sourceW, dx)
    error('Source width sourceW is not a multiple of dx');
end
if rem(sourceL, dx)
    error('Source length sourceL is not a multiple of dx');
end

% Number of grid points to represent element size
Nwidth  = round(sourceW / dx);
Nlength = round(sourceL / dx);

% Create mask in top left corner
mask = zeros(Ny, Nx);
mask(1:Nlength, 1:Nwidth ) = ones(Nlength, Nwidth);

% Calculate the current centre location of the element
ixc = round(sourceW / (2 * dx));
iyc = round(sourceL / (2 * dx)); 

% Calculate the shift needed to move element to the required location
shift_x = ix0 - ixc;
shift_y = iy0 - iyc;

% Move the element to the required location
mask = circshift(mask, [shift_y, shift_x]);

end