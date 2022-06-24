function rgb_bins = getRGBbins(wl_start,wl_end,nBins)
% GETRGBBINS computes the RGB value per bin for plotting 'true' colour
%   
%   Usage:
%   rgb_bins = getRGBbins(wl_start,wl_end,nBins,nDims) computes the
%   N x 3 matrix rgb_bins_m containing the RGB values for N bins. 
%   'wl_start' is the wavelength of the first bin, 'wl_end' of the last bin
%   and the number of bins N is 'nBins'. Wavelengths are in nm.
%
%   SEE ALSO: PLOTTRUECOLOUR

[rgb,lambda] = spectrumColors;
binSizeHalf = 0.5*(wl_end - wl_start)/(nBins - 1);
wl_bins = linspace(wl_start,wl_end,nBins);

% mean over bins in binSize
rgb_bins = nan(nBins,3);
for ii = 1:numel(wl_bins)
    idx = floor(wl_bins(ii)-binSizeHalf:wl_bins(ii)+binSizeHalf); % get idx of bins
    rgb_bins(ii,:) = sum(rgb(ismember(lambda',idx),:),1)./numel(idx); % take mean over bins. Bins outside range lambda are counted as [0 0 0]
end

end