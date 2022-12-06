function mask = createElementMask(Nx, Ny, dx, sourceW, sourceL, ix0, iy0)

% ix0 is an index estimate for the element centre x-cooridnate in dx space
% iy0 is an index estimate for the element centre y-cooridnate in dx space

% builds a 'continuous' mask

% Find greatest common divisor of the input values - the upsampled grid
% spacing (rounded by 1e9 to make everythin integer - removed later)
K = 1e9 * [dx; sourceW / 2; sourceL / 2];
dx_up = K(1);
for i = 2:numel(K)
    dx_up = gcd(dx_up, K(i));
end

% Calculate parameters
dx_up           = dx_up / 1e9;
upsample_factor = dx / dx_up;
Nwidth          = round(sourceW / dx_up);
Nlength         = round(sourceL / dx_up);
Ny_up           = Ny * upsample_factor;
Nx_up           = Nx * upsample_factor;

% find ix0 and iy0 in upsampled space
ix0_up = ix0 * upsample_factor;
iy0_up = iy0 * upsample_factor;

% Calculate grid vectors (for plotting during development)
x_vec_up = 0:dx_up:(Nx_up-1) * dx_up;
x_vec    = x_vec_up(1:upsample_factor:end) + (dx_up * ((upsample_factor-1)/2));

% Generate source mask in upsampled-space from known dimensions
mask_up = zeros(Ny_up, Nx_up);
mask_up(1:Nlength, 1:Nwidth ) = ones(Nlength, Nwidth);

% Calculate the current centre location of the element
ixc = round(sourceW / (2 * dx_up));
iyc = round(sourceL / (2 * dx_up)); 

% Calculate the shift needed to move element to the required location
shift_x = ix0_up - ixc;
shift_y = iy0_up - iyc;

% Move the element to the required location
mask_up = circshift(mask_up, [shift_y, shift_x]);

% Resize to achieve the required grid spacing
mask = imresize(mask_up, [Ny, Nx], 'bilinear', 'Antialiasing', true);

return
% 
% mask_up2  = imresize(mask_down, [Ny_up, Nx_up], 'bilinear');
% 
% diff = mask_up - mask_up2;
% 
% figure;
% subplot(3, 1, 1);
% imagesc(mask_up);
% axis image
% 
% subplot(3, 1, 2);
% imagesc(mask_down);
% axis image
% 
% subplot(3, 1, 3);
% hold on;
% plot(x_vec_up, mask_up(round(Ny_up/2),:), 'k');
% % plot(x_vec_up, mask_up2(round(Ny_up/2),:), 'b');
% plot(x_vec, mask_down(round(Ny/2),:), 'b');
% % downsample again
% 
% test_vec = mask_down(round(Ny/2), :);
% test_vec_up = mask_up(round(Ny_up/2), :);
% shift = 0.05e-3;
% 
% figure;
% for idx = 1:100
%     shift_pts = (idx * shift) / dx;
%     shift_pts_up = round((idx * shift) / dx_up);
% 
%     test_vec    = mask_down(round(Ny/2), :);
%     test_vec_up = mask_up(round(Ny_up/2), :);
% 
%     B = fracCircShift(test_vec, shift_pts);
%     B_up = circshift(test_vec_up, shift_pts_up);
% 
%     clf(gcf);
%     hold on
%     plot(x_vec, B, 'b');
%     plot(x_vec_up, B_up, 'k')
%     drawnow
%     pause(0.2);
% end
% 
% end