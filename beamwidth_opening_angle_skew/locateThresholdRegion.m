function [J, I, K] = locateThresholdRegion(volume_data, thresh, options)
%LOCATETHRESHOLDREGION find the centroid, edges, and width of a thresholded region in a 3D array
%
% DESCRIPTION:
%     locateThresholdRegion is designed for 3D datasets. A threshold is
%     applied to define the region of interest, and the following region
%     properties are calculated: weighted centroid, bounding box edges, and
%     bounding box width.
%
% USAGE:
%     [J, I, K] = locateThresholdRegion(volume, thresh)
%
% INPUTS:
%     volume_data - [numeric] 3D array
%     thresh      - [numeric] threshold in decibels, should be negative
%                   [dB]
%
% OPTIONAL INPUTS:
%     Round       - [boolean] whether to round the ic, imin, imax fields
%                   before returning them
%
% OUTPUS:
%     J           - [structure] data for columns  (dim 2)
%     I           - [structure] data for rows     (dim 1)
%     K           - [structure] data for planes   (dim 3)
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

% Kind the weighted centroid
stats = regionprops3( mask, volume_data, "WeightedCentroid", 'BoundingBox');
I.ic  = stats.WeightedCentroid(2);
J.ic  = stats.WeightedCentroid(1);
K.ic  = stats.WeightedCentroid(3);

% Kind the bounding box edges and widths
I.imin = stats.BoundingBox(2);
I.imax = stats.BoundingBox(2) + stats.BoundingBox(5);
I.w    = stats.BoundingBox(5);
J.imin = stats.BoundingBox(1);
J.imax = stats.BoundingBox(1) + stats.BoundingBox(4);
J.w    = stats.BoundingBox(4);
K.imin = stats.BoundingBox(3);
K.imax = stats.BoundingBox(3) + stats.BoundingBox(6);
K.w    = stats.BoundingBox(6);

if options.Round
    I.ic = round( I.ic );
    J.ic = round( J.ic );
    K.ic = round( K.ic );

    I.imin = round( I.imin );
    J.imin = round( J.imin );
    K.imin = round( K.imin );

    I.imax = round( I.imax );
    J.imax = round( J.imax );
    K.imax = round( K.imax );

end