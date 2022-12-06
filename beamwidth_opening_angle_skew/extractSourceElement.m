function extract_p = extractSourceElement(source_p, dx, i0, options)
%EXTRACTSOURCEELEMENT Extract source pressure for single element from array
%
% DESCRIPTION:
%     extractSourceElement extracts the source pressure distribution for a
%     single element, from a source plane containing the pressure
%     distribution of an entire array. First the pressure squared integral
%     is calculated and fourier upsampled. A binary mask is constructed
%     using the known physical size of the elements. The user provides an
%     initial guess for the centre of the element, and then the actual
%     centre of the element is found using 2D crosscorrelation between the
%     mask and the pressure squared integral. The mask is shifted to the
%     optimal position and then re-sampled to the original grid spacing
%     before being applied to the source pressure in 3D to perform the
%     extraction of a single element from the array.
%
% USAGE:
%     extract_p = extractSourceElement(source_p, dx, i0)
%
% INPUTS:
%     source_p     - [numeric] 3D pressure array with size (Ny, Nx, Nt)
%                    containing the source pressure for the entire array
%                    [Pa]
%     dx           - [numeric] grid spacing of the input plane [m]
%     i0           - [numeric] two-element [row (y), col (x)] array
%                    containing the initial estimate of the element
%                    centroid [pts]
%
% OPTIONAL INPUTS:
%     Border       - [integer] two-element [row (y), col (x)] array
%                    containing the number of border pixels to remove
%                    before upsampling (restored at the end) [pts]
%     UpSample     - [integer] integer value representing the upsample
%                    factor to use (must be odd)
%     ExtraPlot    - [boolean] whether to produce figures showing the
%                    progress of the extraction procedure
%     PhysicalSize - [numeric] two-element [width (x), length (y)] array
%                    containing the physical size of the elements. This is
%                    used for fine adjustment of the mask position.
%                    Alignment performance is poor if these values are
%                    overestimated. [m]
%     MaskSize     - [numeric] two-element [width (x), length (y)] array
%                    containing the physical size of the mask used to
%                    extract the pressure [m]. Should generally be larger
%                    than PhysicalSize.
%
% OUTPUTS:
%     extract_p    - [numeric] 3D pressure array with size (Ny, Nx, Nt)
%                    containing the source pressure for a single element
%                    [Pa]
%
% ABOUT:
%     author       - Morgan Roberts
%     date         - 6/12/22


arguments
    source_p
    dx
    i0
    options.Border       = [10, 80];
    options.UpSample     = 7; % MUST BE ODD
    options.ExtraPlot    = true;
    options.PhysicalSize = [0.9e-3, 9e-3];
    options.MaskSize     = [2.1e-3, 11.1e-3];
end

Nt = size(source_p, 3);

% Create pressure squared integral matrix for mask alignment
ssp = sum(source_p.^2, 3);

% remove border to save memory after upsampling
x_border = options.Border(2);
y_border = options.Border(1);
ssp_cut  = ssp( y_border+1:end-y_border, x_border+1:end-x_border );

% Estimated element centroid location [pts]
ix0 = i0(2) - x_border;
iy0 = i0(1) - y_border;

% Size of extracted matrix
Nx_cut = size(ssp_cut, 2);
Ny_cut = size(ssp_cut, 1);

% Fourier upsample the pressure squared integral to allow detection of the
% element centroid
ups    = options.UpSample;
sz     = size(ssp_cut);
sz_up  = sz * options.UpSample;
ssp_up = interpftn(ssp_cut, sz_up);

% Convert the grid parameters to up-sampled space
ix0_up = ix0 * ups;
iy0_up = iy0 * ups;
dx_up  = dx / ups;
Nx_up  = Nx_cut * ups;
Ny_up  = Ny_cut * ups;

% Create a binary mask in upsampled space using known physical dimensions
% of the elements.
w_phys  = options.PhysicalSize(1);
l_phys  = options.PhysicalSize(2);
mask_up = createElementMask(Nx_up, Ny_up, dx_up, w_phys, l_phys, ix0_up, iy0_up);

% Find the actual centroid of the element in upsampled space using manual
% implementation of 2d cross correlation (in a restricted region) between
% the mask and the pressure squared integral
max_shift = 15; % [pts]
shifts    = -max_shift:max_shift;
Npos      = length(shifts);
r         = zeros(Npos, Npos);
for xdx = 1:Npos
    for ydx = 1:Npos
        x_shift     = shifts(xdx);
        y_shift     = shifts(ydx);
        mask_shift  = circshift(mask_up, [y_shift, x_shift]);
        r(ydx, xdx) = sum(mask_shift .* ssp_up, 'all');
    end
end

% Find the optimal shift to align mask with pressure squared integral
[~,I]          = max(r, [], 'all');
[y_opt, x_opt] = ind2sub(size(r), I);
x_shift        = shifts(x_opt);
y_shift        = shifts(y_opt);

% Create the final mask using the required mask dimensions
w_mask   = options.MaskSize(1);
l_mask   = options.MaskSize(2);
x_opt    = ix0_up + x_shift;
y_opt    = iy0_up + y_shift;
mask_opt = createElementMask(Nx_up, Ny_up, dx_up, w_mask, l_mask, x_opt, y_opt);

% Apply correction (imresize and interpftn handle the physical location of
% upsampled grid points differently)
mask_corr = circshift(mask_opt, [floor(ups / 2), floor(ups / 2)]);

% resample to produce the aligned mask interpolated to the original grid
mask = imresize(mask_corr, sz, 'bilinear', 'Antialiasing', true);

% Restore the mask to the original size of the input matrix
mask_2d = zeros(size(ssp));
mask_2d( y_border+1:end-y_border, x_border+1:end-x_border ) = mask;

% Replicate mask for all time steps
mask_3d = repmat(mask_2d, [1, 1, Nt]);

% Apply mask to source plane to extract the pressure from a single element
extract_p = mask_3d .* source_p;




% If necessary, plot figures to show the process
if options.ExtraPlot

    % Visualise the original and upsampled pressure squared integrals
    figure;
    subplot(3, 1, 1);
    imagesc(ssp_cut);
    axis image
    xlabel('x-position [samples]');
    ylabel('y-position [samples]');
    c = colorbar;
    ylabel(c, 'Pressure Squared Integral');
    colormap(getBatlow)
    title('Pressure squared integral');
    
    subplot(3, 1, 2);
    imagesc(ssp_up);
    axis image
    xlabel('x-position [samples]');
    ylabel('y-position [samples]');
    c = colorbar;
    ylabel(c, 'Pressure Squared Integral');
    colormap(getBatlow)
    title('Upsampled pressure squared integral');
    
    subplot(3, 1, 3);
    hold on;
    plot( 1:ups:ups*size(ssp_cut, 2), ssp_cut(iy0,:), 'kx' );
    plot( ssp_up( ( (iy0 - 1) * ups) + 1,: ), 'b' )
    xlabel('x-position [samples]');
    ylabel('Pressure Squared Integral');
    legend({'Original Data', 'Upsampled Data'})

% ------------------------------------------------------------------------

    % Visualise the mask generation
    % Re-sample interpolate the mask back to the original grid size
    mask_check = imresize(mask_up, sz, 'bilinear', 'Antialiasing', true);
    
    figure;
    subplot(3, 1, 1);
    imagesc(mask_up);
    axis image
    xlabel('x-position [samples]');
    ylabel('y-position [samples]');
    title('Upsampled mask (binary)');
    
    subplot(3, 1, 2);
    imagesc(mask_check);
    axis image
    xlabel('x-position [samples]');
    ylabel('y-position [samples]');
    title('Re-sampled mask (interpolated)');
    
    i_vec = (1:ups:Nx_up) + ((ups - 1) / 2);
    
    subplot(3, 1, 3);
    hold on;
    plot(i_vec, mask_check(iy0,:), 'k');
    plot(mask_up(iy0_up,: ), 'b');
    xlabel('x-position [samples]');
    ylabel('Mask Value');
    legend({'Re-sampled Mask', 'Upsampled Mask'})

% ------------------------------------------------------------------------    

    % Plot the mask alignment process 
    figure;
    imagesc(shifts, shifts, r);
    c = colorbar;
    axis image
    colormap(getBatlow);
    hold on;
    plot(x_shift, y_shift, 'kx');
    xlabel('x-shift [pts]');
    ylabel('y-shift [pts]');
    ylabel(c, 'Cross correlation value');
    title('Mask alignment ')


% ------------------------------------------------------------------------  
    
    % Plot the masked pressure
    figure;

    subplot(2, 2, 1);
    hold on;
    plot(mask_opt(iy0_up,:) * max(ssp_up( ( (iy0 - 1) * ups) + 1,: )), 'k');
    plot(ssp_up( ( (iy0 - 1) * ups) + 1,: ), 'b')
    xlabel('x-position [pts]');
    ylabel('Pressure squared integral');
    legend({'Mask', 'Pressure squared integral'});
    title('Upsampled Masking');

    subplot(2, 2, 2);
    hold on;
    plot(mask(iy0,:) * max(ssp_cut(iy0,:)), 'k');
    plot(ssp_cut(iy0,:), 'b');
    xlabel('x-position [pts]');
    ylabel('Pressure squared integral');
    legend({'Mask', 'Pressure squared integral'});    
    title('Resampled Masking');    
    
    subplot(2, 2, 3);
    imagesc((mask_opt + 0.20) .* ssp_up);
    axis image
    colormap(getBatlow);   
    c = colorbar;
    xlabel('x-position [pts]');   
    ylabel('y-position [pts]');       
    ylabel(c, 'Pressure Squared Integral');
    
    subplot(2, 2, 4);
    imagesc((mask + 0.20) .* ssp_cut);
    axis image
    colormap(getBatlow);
    c = colorbar;
    xlabel('x-position [pts]');   
    ylabel('y-position [pts]');       
    ylabel(c, 'Pressure Squared Integral');    


end

end