function plotTrueColor2D(obj,ha,id,Zidx,Tidx,varargin)
%PLOTTRUECOLOR plots a 'true' color image
%   Wrapper for plotRGBbins
% 
%   Usage:
%   plotTrueColor2D(obj,ha,id,Zidx,Tidx,rgb_bins) plot a 'true' color image
%   in axis handle ha (leave empty for current open axis), for image with
%   index id, and Zidx for the index in the Z-direction and Tidx for the
%   index in the T-direction. rgb_bins is the output of function
%   getRGBbins.
%   plotTrueColor2D(obj,ha,id,Zidx,Tidx,wl_start,wl_end,nBins) is same as
%   above, rgb_bins is computed based on the parameters wl_start,wl_end,
%   and nBins.
%   'wl_start' is the wavelength of the first bin, 'wl_end' of the last bin
%   and the number of bins N is 'nBins'. Wavelengths are in nm.

if isempty(ha)
    ha = gca;
end

if isempty(id)
    id = 1;
end

if isempty(Zidx)
    Zidx = 1;
end

if isempty(Tidx)
    Tidx = 1;
end

I = obj.getImageProcessed(id,'reshaped','x','y','z','t','c');
I = squeeze(I(:,:,Zidx,Tidx,:)); % [X Y C]

if nargin == 6
    rgb_bins = varargin{1};
    plotTrueColour(ha,I,rgb_bins)
elseif nargin == 9
    wl_start = varargin{1};
    wl_end = varargin{2};
    nBins = varargin{3};
    plotTrueColour(ha,I,getRGBbins(wl_start,wl_end,nBins));
else
    error('Number of inputs not supported. Got %i, need 6 or 9.',nargin)
end

