function [X, Y, F] = locateThresholdRegion(volume_data, thresh, options)
%LOCATETHRESHOLDREGION find the centroid, edges, and width of a thresholded region in a 3D array
%
% DESCRIPTION:
%     locateThresholdRegion is designed for 3D datasets containing
%     amplitude pressure spectra over a XY plane. A threshold is applied to
%     define the region of interest, and the following region properties
%     are calculated:weighted centroid, bounding box edges, and bounding
%     box width.
%
% USAGE:
%     [X, Y, F] = locateThresholdRegion(volume, thresh)
%
% INPUTS:
%     volume_data - [numeric] pressure array with size (Ny, Nx, Nf)
%                   containing amplitude spectra for each point over a XY
%                   plane [Pa]
%     thresh      - [numeric] threshold in decibels, should be negative
%                   [dB]
%
% OPTIONAL INPUTS:
%     Round       - [boolean] whether to round the ic, imin, imax fields
%                   before returning them
%
% OUTPUS:
%     X           - [structure] data for dimension 2 (x-direction)
%     Y           - [structure] data for dimension 1 (y-direction)
%     F           - [structure] data for dimension 3 (f-direction)
%
% FIELDS:
%     ic          - index of centroid
%     imin        - index of bounding box start
%     imin        - index of bounding box end
%     w           - width of bounding box
%
% ABOUT:
%     author      - Morgan Roberts
%     date        - 7/12/22

arguments
    volume_data
    thresh
    options.Round = true;
end

% Convert to dB and create mask above threshold
volume_dB = 20 * log10( volume_data / max( volume_data(:) ) );
mask      = volume_dB > thresh;

% Keep only the largest connected component
CC             = bwconncomp(mask);
numOfPixels    = cellfun(@numel,CC.PixelIdxList);
[~,indexOfMax] = max(numOfPixels);
mask           = zeros(size(mask));
mask(CC.PixelIdxList{indexOfMax}) = 1;

% Find the weighted centroid
stats = regionprops3( mask, volume_data, "WeightedCentroid", 'BoundingBox');
Y.ic  = stats.WeightedCentroid(2);
X.ic  = stats.WeightedCentroid(1);
F.ic  = stats.WeightedCentroid(3);

% Find the bounding box edges and widths
Y.imin = stats.BoundingBox(2);
Y.imax = stats.BoundingBox(2) + stats.BoundingBox(5);
Y.w    = stats.BoundingBox(5);
X.imin = stats.BoundingBox(1);
X.imax = stats.BoundingBox(1) + stats.BoundingBox(4);
X.w    = stats.BoundingBox(4);
F.imin = stats.BoundingBox(3);
F.imax = stats.BoundingBox(3) + stats.BoundingBox(6);
F.w    = stats.BoundingBox(6);

if options.Round
    Y.ic = round( Y.ic );
    X.ic = round( X.ic );
    F.ic = round( F.ic );

    Y.imin = round( Y.imin );
    X.imin = round( X.imin );
    F.imin = round( F.imin );

    Y.imax = round( Y.imax );
    X.imax = round( X.imax );
    F.imax = round( F.imax );

end