function plotTrueColour(ha,I,rgb_bins)
%PLOTRUECOLOR Plots image with true colours
%   
%   Usage:
%   plotTrueColour(I,rgb_bins) plots 2-D or 3-D image I, which has
%   'c' as last dimension using the true colour values in rgb_bins.
%   Use function getRGBbins to compute rgb_bins.
%   Use I = getImageProcessed(cl,ii,'mean','x','y','c'); to get image
%
%   SEE ALSO: GETRGBBINS, GETIMAGEPROCESSED

if isempty(ha)
    ha = gca;
end

sz = size(I);
if sz(end) ~= size(rgb_bins,1)
    error('The number of bins in ''rgb_bins'' and I do not match. Found %i RGB bins and %i bins in image.',size(rgb_bins,1),sz(end))
end

nDims = numel(sz);

% permute so that last dimension (here 3rd) is multiplied with rgb values
% per channel
idx = 1:nDims+1;
idx = circshift(idx,nDims-1);
rgb_bins_m = permute(rgb_bins,idx);

% apply rgb to image
I_rgb = squeeze(sum(double(I).*rgb_bins_m,3));
I_rgb = I_rgb./max(I_rgb,[],'all');

imshow(I_rgb,'Parent',ha)
axis image

end

